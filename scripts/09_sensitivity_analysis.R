# ==============================================================================
# Script: 09_sensitivity_analysis.R
# Purpose: Calculate sensitivity/specificity for each significant DMR gene
#          for pancreatic cancer prediction
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(bsseq)
  library(data.table)
  library(pROC)
  library(ggplot2)
  library(GenomicRanges)
  library(httr)
  library(jsonlite)
})

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# Configuration
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
DMRSEQ_DIR <- file.path(paths$output_dir, "dss_results", "dmrseq", "20251228_cut0.05_minR3_bp1000")

logit("=== Sensitivity Analysis for Pancreatic Cancer Prediction ===")

# Load BSseq object
bs_file <- file.path(H5_DIR, "bsseq.rds")
targets_file <- file.path(H5_DIR, "targets.rds")

BSobj <- readRDS(bs_file)
targets <- readRDS(targets_file)

# Load significant DMRs
sig_dmr_file <- file.path(DMRSEQ_DIR, "dmrseq_Significant_DMRs.tsv")
sig_dmrs <- fread(sig_dmr_file)

logit("Loaded ", nrow(sig_dmrs), " significant DMRs")

# Ensembl REST API for gene annotation
get_gene_at_position <- function(chr, start, end) {
  chr_clean <- gsub("chr", "", chr)
  url <- paste0("https://rest.ensembl.org/overlap/region/human/", 
                chr_clean, ":", start, "-", end, 
                "?feature=gene;content-type=application/json")
  tryCatch({
    response <- GET(url, timeout(10))
    if(status_code(response) == 200) {
      genes <- fromJSON(content(response, "text", encoding = "UTF-8"))
      if(length(genes) > 0 && "external_name" %in% names(genes)) {
        valid_names <- genes$external_name[!is.na(genes$external_name) & genes$external_name != ""]
        if(length(valid_names) > 0) return(paste(unique(valid_names), collapse = "/"))
      }
    }
    return(NA)
  }, error = function(e) return(NA))
}

# Create GRanges for DMRs
dmr_gr <- GRanges(seqnames = sig_dmrs$seqnames, 
                  ranges = IRanges(start = sig_dmrs$start, end = sig_dmrs$end))
bs_gr <- granges(BSobj)

# Calculate methylation for each DMR and perform ROC analysis
logit("Calculating methylation levels and ROC analysis...")

results <- list()

for(i in 1:nrow(sig_dmrs)) {
  # Get overlapping CpGs
  overlaps <- findOverlaps(dmr_gr[i], bs_gr)
  cpg_idx <- subjectHits(overlaps)
  
  if(length(cpg_idx) > 0) {
    # Get methylation values
    meth_vals <- as.matrix(getMeth(BSobj[cpg_idx, ], type = "raw"))
    mean_meth <- colMeans(meth_vals, na.rm = TRUE)
    
    # Create data frame for ROC
    roc_data <- data.frame(
      methylation = mean_meth,
      group = targets$Group  # "Case" = Cancer, "Control" = Normal
    )
    roc_data <- roc_data[!is.na(roc_data$methylation), ]
    
    if(nrow(roc_data) > 10) {
      # Perform ROC analysis
      # For hypermethylation (beta > 0): higher methylation = cancer
      # For hypomethylation (beta < 0): lower methylation = cancer
      direction <- ifelse(sig_dmrs$beta[i] > 0, "greater", "less")
      
      roc_obj <- tryCatch({
        roc(roc_data$group, roc_data$methylation, 
            levels = c("Control", "Case"), 
            direction = ifelse(direction == "greater", "<", ">"))
      }, error = function(e) NULL)
      
      if(!is.null(roc_obj)) {
        # Find optimal threshold (Youden's J)
        coords_best <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
        
        # Get gene name
        gene_name <- get_gene_at_position(sig_dmrs$seqnames[i], sig_dmrs$start[i], sig_dmrs$end[i])
        if(is.na(gene_name)) gene_name <- paste0(sig_dmrs$seqnames[i], ":", sig_dmrs$start[i])
        
        results[[i]] <- data.frame(
          Rank = i,
          Gene = gene_name,
          Chr = sig_dmrs$seqnames[i],
          Start = sig_dmrs$start[i],
          End = sig_dmrs$end[i],
          Direction = ifelse(sig_dmrs$beta[i] > 0, "Hypermethylation", "Hypomethylation"),
          Beta = round(sig_dmrs$beta[i], 4),
          Qvalue = signif(sig_dmrs$qval[i], 3),
          AUC = round(auc(roc_obj), 4),
          Sensitivity = round(coords_best$sensitivity, 4),
          Specificity = round(coords_best$specificity, 4),
          Threshold = round(coords_best$threshold, 4),
          PPV = NA,  # Will calculate below
          NPV = NA
        )
        
        # Calculate PPV and NPV
        n_case <- sum(roc_data$group == "Case")
        n_control <- sum(roc_data$group == "Control")
        prevalence <- n_case / (n_case + n_control)
        
        sens <- coords_best$sensitivity
        spec <- coords_best$specificity
        
        # PPV = (Sens * Prev) / (Sens * Prev + (1-Spec) * (1-Prev))
        ppv <- (sens * prevalence) / (sens * prevalence + (1 - spec) * (1 - prevalence))
        # NPV = (Spec * (1-Prev)) / (Spec * (1-Prev) + (1-Sens) * Prev)
        npv <- (spec * (1 - prevalence)) / (spec * (1 - prevalence) + (1 - sens) * prevalence)
        
        results[[i]]$PPV <- round(ppv, 4)
        results[[i]]$NPV <- round(npv, 4)
        
        Sys.sleep(0.1)  # API rate limit
      }
    }
  }
}

# Combine results
results_df <- do.call(rbind, results[!sapply(results, is.null)])

# Sort by Sensitivity (descending)
results_df <- results_df[order(-results_df$Sensitivity), ]
results_df$Sensitivity_Rank <- 1:nrow(results_df)

# Reorder columns
results_df <- results_df[, c("Sensitivity_Rank", "Gene", "Direction", "Sensitivity", "Specificity", 
                              "AUC", "PPV", "NPV", "Threshold", "Beta", "Qvalue", "Chr", "Start", "End")]

logit("Analysis complete for ", nrow(results_df), " DMRs")

# Save results
out_file <- file.path(DMRSEQ_DIR, "sensitivity_analysis.tsv")
fwrite(results_df, out_file, sep = "\t")
logit("Results saved to: ", out_file)

# Print top results
logit("")
logit("=== Top DMRs by Sensitivity for Pancreatic Cancer Prediction ===")
logit("")

cat(sprintf("%-4s %-30s %-15s %-10s %-10s %-8s\n", 
            "Rank", "Gene", "Direction", "Sensitivity", "Specificity", "AUC"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for(i in 1:min(nrow(results_df), 33)) {
  cat(sprintf("%-4d %-30s %-15s %-10.1f%% %-10.1f%% %-8.3f\n",
              results_df$Sensitivity_Rank[i],
              substr(results_df$Gene[i], 1, 30),
              results_df$Direction[i],
              results_df$Sensitivity[i] * 100,
              results_df$Specificity[i] * 100,
              results_df$AUC[i]))
}

logit("")
logit("=== Interpretation ===")
logit("Sensitivity: % of cancer patients correctly identified by this methylation marker")
logit("Specificity: % of normal individuals correctly identified")
logit("AUC: Overall discriminative ability (0.5 = random, 1.0 = perfect)")
logit("")

# Create visualization
logit("Creating sensitivity plot...")

plot_df <- results_df[1:min(nrow(results_df), 20), ]
# Handle duplicate gene names by adding rank
plot_df$Gene_Label <- paste0(plot_df$Sensitivity_Rank, ". ", plot_df$Gene)
plot_df$Gene_Label <- factor(plot_df$Gene_Label, levels = rev(plot_df$Gene_Label))

p <- ggplot(plot_df, aes(x = Gene_Label, y = Sensitivity * 100, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Sensitivity * 100, 1), "%")), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Hypermethylation" = "#D73027", "Hypomethylation" = "#4575B4")) +
  labs(title = "Pancreatic Cancer Prediction Sensitivity by DMR Gene",
       subtitle = paste0("Top ", nrow(plot_df), " genes ranked by sensitivity (n=94 patients, n=50 normals)"),
       x = "Gene",
       y = "Sensitivity (%)",
       fill = "Methylation Direction") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(0, 110)

ggsave(file.path(DMRSEQ_DIR, "plots", "07_Sensitivity_Ranking.png"), p, 
       width = 12, height = 10, dpi = 300)

logit("Plot saved to: ", file.path(DMRSEQ_DIR, "plots", "07_Sensitivity_Ranking.png"))
logit("=== Analysis Complete ===")
