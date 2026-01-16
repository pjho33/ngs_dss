# ==============================================================================
# Script: 07b_dmrseq_visualization.R
# Purpose: Visualize dmrseq results (Heatmap, PCA, Volcano, etc.)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(bsseq)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(GenomicRanges)
  library(httr)
  library(jsonlite)
})

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# Configuration
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
DMRSEQ_BASE <- file.path(paths$output_dir, "dss_results", "dmrseq")

# BSseq 객체 로드
bs_file <- file.path(H5_DIR, "bsseq.rds")
targets_file <- file.path(H5_DIR, "targets.rds")

if(file.exists(bs_file)) {
  BSobj <- readRDS(bs_file)
  targets <- readRDS(targets_file)
} else {
  stop("BSseq object not found!")
}

# Ensembl REST API로 유전자 조회 함수
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

# 모든 dmrseq 결과 폴더 찾기
subdirs <- list.dirs(DMRSEQ_BASE, recursive = FALSE, full.names = TRUE)

logit("Found ", length(subdirs), " dmrseq result folders")

for(subdir in subdirs) {
  logit("========================================")
  logit("Processing: ", basename(subdir))
  
  dmr_file <- file.path(subdir, "dmrseq_DMRs.tsv")
  sig_dmr_file <- file.path(subdir, "dmrseq_Significant_DMRs.tsv")
  plot_dir <- file.path(subdir, "plots")
  
  if(!file.exists(dmr_file)) {
    logit("  -> Skipping: No DMR file found")
    next
  }
  
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  dmrs <- fread(dmr_file)
  sig_dmrs <- if(file.exists(sig_dmr_file)) fread(sig_dmr_file) else dmrs[qval < 0.05]
  
  logit("  Loaded: ", nrow(dmrs), " total DMRs, ", nrow(sig_dmrs), " significant")
  
  # ==============================================================================
  # 1. Volcano Plot
  # ==============================================================================
  logit("  1. Volcano Plot...")
  
  dmrs[, neg_log_q := -log10(qval)]
  dmrs[, neg_log_q := pmin(neg_log_q, 50)]
  dmrs[, significance := "Not Significant"]
  dmrs[qval < 0.05 & beta > 0, significance := "Hyper (Case)"]
  dmrs[qval < 0.05 & beta < 0, significance := "Hypo (Case)"]
  
  p_volcano <- ggplot(dmrs, aes(x = beta, y = neg_log_q, color = significance)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Hyper (Case)" = "red", "Hypo (Case)" = "blue", "Not Significant" = "grey60")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "dmrseq: Volcano Plot",
         x = "Methylation Difference (beta: Case - Control)",
         y = "-log10(q-value)",
         color = "Status") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, "01_Volcano_Plot.png"), p_volcano, width = 10, height = 8, dpi = 300)
  
  # ==============================================================================
  # 2. Hyper/Hypo Pie Chart
  # ==============================================================================
  logit("  2. Pie Chart...")
  
  dmrs[, direction := ifelse(beta > 0, "Hypermethylated", "Hypomethylated")]
  direction_counts <- dmrs[, .N, by = direction]
  direction_counts[, pct := round(N / sum(N) * 100, 1)]
  direction_counts[, label := paste0(direction, "\n(n=", N, ", ", pct, "%)")]
  
  p_pie <- ggplot(direction_counts, aes(x = "", y = N, fill = direction)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = c("Hypermethylated" = "#E41A1C", "Hypomethylated" = "#377EB8")) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = "dmrseq: DMR Direction (Case vs Control)") +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(file.path(plot_dir, "02_Direction_PieChart.png"), p_pie, width = 8, height = 8, dpi = 300)
  
  # ==============================================================================
  # 3. DMR Heatmap (Significant DMRs)
  # ==============================================================================
  if(nrow(sig_dmrs) > 0) {
    logit("  3. DMR Heatmap...")
    
    dmr_gr <- GRanges(seqnames = sig_dmrs$seqnames, 
                      ranges = IRanges(start = sig_dmrs$start, end = sig_dmrs$end))
    bs_gr <- granges(BSobj)
    
    dmr_meth_list <- list()
    for(i in 1:length(dmr_gr)) {
      overlaps <- findOverlaps(dmr_gr[i], bs_gr)
      cpg_idx <- subjectHits(overlaps)
      if(length(cpg_idx) > 0) {
        meth_vals <- as.matrix(getMeth(BSobj[cpg_idx, ], type = "raw"))
        dmr_meth_list[[i]] <- colMeans(meth_vals, na.rm = TRUE)
      } else {
        dmr_meth_list[[i]] <- rep(NA, ncol(BSobj))
      }
    }
    
    dmr_meth_mat <- do.call(rbind, dmr_meth_list)
    
    # Gene annotation
    logit("     -> Annotating DMRs with gene names...")
    gene_names <- character(nrow(sig_dmrs))
    for(i in 1:nrow(sig_dmrs)) {
      gene_name <- get_gene_at_position(sig_dmrs$seqnames[i], sig_dmrs$start[i], sig_dmrs$end[i])
      gene_names[i] <- if(!is.na(gene_name) && gene_name != "") gene_name else paste0(sig_dmrs$seqnames[i], ":", sig_dmrs$start[i])
      Sys.sleep(0.1)
    }
    
    if(any(duplicated(gene_names))) {
      dup_idx <- which(duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE))
      for(idx in dup_idx) {
        gene_names[idx] <- paste0(gene_names[idx], " (", sig_dmrs$seqnames[idx], ":", sig_dmrs$start[idx], ")")
      }
    }
    
    rownames(dmr_meth_mat) <- gene_names
    colnames(dmr_meth_mat) <- colnames(BSobj)
    
    # NA 처리
    for(i in 1:nrow(dmr_meth_mat)) {
      na_idx <- is.na(dmr_meth_mat[i, ])
      if(any(na_idx)) {
        row_mean <- mean(dmr_meth_mat[i, !na_idx], na.rm = TRUE)
        if(is.na(row_mean)) row_mean <- 0.5
        dmr_meth_mat[i, na_idx] <- row_mean
      }
    }
    dmr_meth_mat[is.na(dmr_meth_mat)] <- 0.5
    
    # 유효한 행만 유지
    valid_rows <- rowSums(is.finite(dmr_meth_mat)) == ncol(dmr_meth_mat)
    dmr_meth_mat <- dmr_meth_mat[valid_rows, , drop = FALSE]
    sig_dmrs_valid <- sig_dmrs[valid_rows, ]
    
    if(nrow(dmr_meth_mat) > 1) {
      sample_order <- order(targets$Group, decreasing = FALSE)
      dmr_meth_ordered <- dmr_meth_mat[, sample_order, drop = FALSE]
      targets_ordered <- targets[sample_order, ]
      
      anno_col <- data.frame(Group = targets_ordered$Group, row.names = targets_ordered$SampleID)
      anno_colors <- list(Group = c("Control" = "#4DAF4A", "Case" = "#E41A1C"))
      
      n_case <- sum(targets_ordered$Group == "Case")
      
      viridis_colors <- colorRampPalette(c("#440154", "#482878", "#3E4A89", "#31688E", 
                                            "#26828E", "#1F9E89", "#35B779", "#6DCD59", 
                                            "#B4DE2C", "#FDE725"))(256)
      
      heatmap_height <- max(600, min(1500, 100 + nrow(dmr_meth_ordered) * 25))
      
      # 3a: Standard heatmap with clustering
      row_vars <- apply(dmr_meth_ordered, 1, var, na.rm = TRUE)
      do_cluster_rows <- all(row_vars > 1e-10, na.rm = TRUE) && nrow(dmr_meth_ordered) > 2
      
      png(file.path(plot_dir, "03a_Heatmap_DMR.png"), 
          width = 1400, height = heatmap_height, res = 150)
      pheatmap(dmr_meth_ordered,
               scale = "none",
               clustering_distance_rows = if(do_cluster_rows) "euclidean" else NA,
               clustering_method = "ward.D2",
               cluster_rows = do_cluster_rows,
               cluster_cols = FALSE,
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               color = viridis_colors,
               main = paste("dmrseq: Significant DMRs (n=", nrow(dmr_meth_ordered), ", qval<0.05)"),
               fontsize_row = 8,
               fontsize_col = 6,
               show_rownames = TRUE,
               show_colnames = FALSE,
               gaps_col = n_case,
               border_color = NA)
      dev.off()
      
      # 3b: Grouped heatmap (Hyper/Hypo separated)
      # dmrseq에서 beta > 0: Case > Control → Hypermethylation
      dmr_direction <- ifelse(sig_dmrs_valid$beta > 0, "Hypermethylation", "Hypomethylation")
      
      hyper_idx <- which(dmr_direction == "Hypermethylation")
      hypo_idx <- which(dmr_direction == "Hypomethylation")
      hyper_order <- hyper_idx[order(-abs(sig_dmrs_valid$beta[hyper_idx]))]
      hypo_order <- hypo_idx[order(-abs(sig_dmrs_valid$beta[hypo_idx]))]
      row_order <- c(hyper_order, hypo_order)
      
      if(length(row_order) > 0) {
        dmr_meth_grouped <- dmr_meth_ordered[row_order, , drop = FALSE]
        dmr_direction_ordered <- dmr_direction[row_order]
        
        anno_row <- data.frame(
          Direction = factor(dmr_direction_ordered, levels = c("Hypermethylation", "Hypomethylation")),
          row.names = rownames(dmr_meth_grouped)
        )
        
        anno_colors_grouped <- list(
          Group = c("Control" = "#4DAF4A", "Case" = "#E41A1C"),
          Direction = c("Hypermethylation" = "#D73027", "Hypomethylation" = "#4575B4")
        )
        
        n_hyper <- length(hyper_order)
        
        png(file.path(plot_dir, "03b_Heatmap_DMR_Grouped.png"), 
            width = 1600, height = heatmap_height, res = 150)
        pheatmap(dmr_meth_grouped,
                 scale = "none",
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 annotation_col = anno_col,
                 annotation_row = anno_row,
                 annotation_colors = anno_colors_grouped,
                 color = viridis_colors,
                 main = paste0("dmrseq: DMRs Grouped by Direction\n",
                              "(Hyper: n=", n_hyper, ", Hypo: n=", length(hypo_order), ")"),
                 fontsize_row = 8,
                 fontsize_col = 6,
                 show_rownames = TRUE,
                 show_colnames = FALSE,
                 gaps_col = n_case,
                 gaps_row = n_hyper,
                 border_color = NA)
        dev.off()
        
        logit("     -> Heatmaps created: ", n_hyper, " Hyper + ", length(hypo_order), " Hypo DMRs")
        
        # ==============================================================================
        # 4. PCA with DMR regions only
        # ==============================================================================
        logit("  4. PCA with DMR regions...")
        
        pca_res_dmr <- prcomp(t(dmr_meth_ordered), scale. = TRUE)
        pca_df_dmr <- data.frame(
          PC1 = pca_res_dmr$x[, 1],
          PC2 = pca_res_dmr$x[, 2],
          Sample = rownames(pca_res_dmr$x),
          Group = targets_ordered$Group
        )
        
        var_explained_dmr <- round(100 * pca_res_dmr$sdev^2 / sum(pca_res_dmr$sdev^2), 1)
        
        p_pca_dmr <- ggplot(pca_df_dmr, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
          geom_point(size = 4) +
          geom_text(vjust = -1, size = 3) +
          scale_color_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
          labs(title = paste0("dmrseq PCA: ", nrow(dmr_meth_ordered), " Significant DMRs"),
               subtitle = paste0("(Hyper: ", n_hyper, ", Hypo: ", length(hypo_order), ")"),
               x = paste0("PC1 (", var_explained_dmr[1], "%)"),
               y = paste0("PC2 (", var_explained_dmr[2], "%)")) +
          theme_bw() +
          theme(legend.position = "bottom")
        
        ggsave(file.path(plot_dir, "04_PCA_DMRs.png"), p_pca_dmr, width = 10, height = 8, dpi = 300)
        logit("     -> PCA created")
      }
    }
  }
  
  # ==============================================================================
  # 5. Chromosome Distribution
  # ==============================================================================
  logit("  5. Chromosome Distribution...")
  
  chr_counts <- dmrs[, .N, by = seqnames]
  chr_counts[, chr_num := gsub("chr", "", seqnames)]
  chr_counts[, chr_num := factor(chr_num, levels = c(1:22, "X", "Y"))]
  chr_counts <- chr_counts[!is.na(chr_num)]
  
  p_chr <- ggplot(chr_counts, aes(x = chr_num, y = N, fill = N)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#FEE08B", high = "#D53E4F") +
    labs(title = "dmrseq: DMR Distribution by Chromosome",
         x = "Chromosome",
         y = "Number of DMRs") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(plot_dir, "05_Chromosome_Distribution.png"), p_chr, width = 10, height = 6, dpi = 300)
  
  # ==============================================================================
  # 6. DMR Width Distribution
  # ==============================================================================
  logit("  6. DMR Width Distribution...")
  
  p_width <- ggplot(dmrs, aes(x = width)) +
    geom_histogram(bins = 50, fill = "#3288BD", color = "white", alpha = 0.8) +
    scale_x_log10() +
    labs(title = "dmrseq: DMR Width Distribution",
         x = "DMR Width (bp, log scale)",
         y = "Count") +
    theme_bw()
  
  ggsave(file.path(plot_dir, "06_Width_Distribution.png"), p_width, width = 8, height = 6, dpi = 300)
  
  logit("  -> Done!")
}

logit("========================================")
logit("All dmrseq visualizations complete!")
