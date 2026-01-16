# ==============================================================================
# Script: 08_visualization_all.R
# Purpose: Run FULL visualization for ALL result subfolders
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(DSS)
  library(bsseq)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(rlang)
  library(scales)
  library(httr)
  library(jsonlite)
})

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# Configuration
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
BASE_OUT_DIR <- file.path(paths$output_dir, "dss_results")

# BSseq 객체 로드
bs_file <- file.path(H5_DIR, "bsseq.rds")
targets_file <- file.path(H5_DIR, "targets.rds")

if(file.exists(bs_file)) {
  BSobj <- readRDS(bs_file)
  targets <- readRDS(targets_file)
} else {
  stop("BSseq object not found!")
}

# 모든 서브폴더 찾기
subdirs <- list.dirs(BASE_OUT_DIR, recursive = FALSE, full.names = TRUE)
subdirs <- subdirs[grepl("^p[0-9]", basename(subdirs))]

logit("Found ", length(subdirs), " result folders to process")

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

for(subdir in subdirs) {
  logit("========================================")
  logit("Processing: ", basename(subdir))
  
  dml_file <- file.path(subdir, "All_DMLs.tsv.gz")
  dmr_file <- file.path(subdir, "Final_DMRs.tsv")
  plot_dir <- file.path(subdir, "plots")
  
  if(!file.exists(dml_file)) {
    logit("  -> Skipping: No DML file found")
    next
  }
  
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  dml_res <- fread(dml_file)
  dmrs <- if(file.exists(dmr_file)) fread(dmr_file) else NULL
  
  logit("  Loaded: ", nrow(dml_res), " DMLs, ", if(!is.null(dmrs)) nrow(dmrs) else 0, " DMRs")
  
  # Parameters from folder name
  folder_params <- strsplit(basename(subdir), "_")[[1]]
  DSS_P_VAL <- as.numeric(gsub("p", "", folder_params[1]))
  DSS_DELTA <- as.numeric(gsub("d", "", folder_params[2]))
  
  # ==============================================================================
  # 1. Volcano Plot
  # ==============================================================================
  logit("  1. Volcano Plot...")
  
  # DSS에서 diff = meanMethy1(Control) - meanMethy2(Case)
  # diff < 0: Case > Control → Hyper (Case)
  # diff > 0: Case < Control → Hypo (Case)
  dml_plot <- copy(dml_res)
  dml_plot[, neg_log_p := -log10(pval)]
  dml_plot[, neg_log_p := pmin(neg_log_p, 50)]
  dml_plot[, significance := "Not Significant"]
  dml_plot[pval < DSS_P_VAL & diff < -DSS_DELTA, significance := "Hyper (Case)"]
  dml_plot[pval < DSS_P_VAL & diff > DSS_DELTA, significance := "Hypo (Case)"]
  
  p_volcano <- ggplot(dml_plot, aes(x = diff, y = neg_log_p, color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Hyper (Case)" = "red", "Hypo (Case)" = "blue", "Not Significant" = "grey60")) +
    geom_hline(yintercept = -log10(DSS_P_VAL), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-DSS_DELTA, DSS_DELTA), linetype = "dashed", color = "black") +
    labs(title = "Volcano Plot: Differentially Methylated Loci",
         x = "Methylation Difference (Control - Case)",
         y = "-log10(p-value)",
         color = "Status") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, "01_Volcano_Plot.png"), p_volcano, width = 10, height = 8, dpi = 300)
  
  # ==============================================================================
  # 2. Manhattan Plot
  # ==============================================================================
  logit("  2. Manhattan Plot...")
  
  chr_order <- c(paste0("chr", c(1:22, "X", "Y")), as.character(c(1:22, "X", "Y")))
  dml_manhattan <- copy(dml_res)
  dml_manhattan <- dml_manhattan[chr %in% chr_order]
  dml_manhattan[, chr_num := gsub("chr", "", chr)]
  dml_manhattan[, chr_num := factor(chr_num, levels = c(1:22, "X", "Y"))]
  
  dml_manhattan <- dml_manhattan[order(chr_num, pos)]
  chr_lengths <- dml_manhattan[, .(max_pos = max(pos)), by = chr_num]
  chr_lengths[, cumsum_pos := cumsum(as.numeric(max_pos)) - max_pos]
  dml_manhattan <- merge(dml_manhattan, chr_lengths[, .(chr_num, cumsum_pos)], by = "chr_num")
  dml_manhattan[, bp_cum := pos + cumsum_pos]
  
  axis_df <- dml_manhattan[, .(center = mean(bp_cum)), by = chr_num]
  dml_manhattan[, chr_color := as.numeric(chr_num) %% 2]
  
  p_manhattan <- ggplot(dml_manhattan, aes(x = bp_cum, y = -log10(pval), color = factor(chr_color))) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c("0" = "#1f78b4", "1" = "#a6cee3"), guide = "none") +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_num) +
    geom_hline(yintercept = -log10(DSS_P_VAL), linetype = "dashed", color = "red") +
    labs(title = "Manhattan Plot: Genome-wide DML Distribution",
         x = "Chromosome", y = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(plot_dir, "02_Manhattan_Plot.png"), p_manhattan, width = 14, height = 6, dpi = 300)
  
  # ==============================================================================
  # 3. PCA Plot
  # ==============================================================================
  logit("  3. PCA Plot...")
  
  set.seed(42)
  n_sites <- min(50000, length(BSobj))
  sample_idx <- sort(sample(1:length(BSobj), n_sites))
  
  M_mat <- as.matrix(getMeth(BSobj[sample_idx, ], type = "raw"))
  M_mat[is.na(M_mat)] <- 0.5
  
  row_vars <- apply(M_mat, 1, var, na.rm = TRUE)
  M_mat <- M_mat[row_vars > 0 & !is.na(row_vars), ]
  
  pca_res <- prcomp(t(M_mat), scale. = TRUE)
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Sample = colnames(M_mat),
    Group = targets$Group
  )
  
  var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3) +
    scale_color_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
    labs(title = "PCA: Sample Clustering by Methylation",
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, "03_PCA_Plot.png"), p_pca, width = 10, height = 8, dpi = 300)
  
  # ==============================================================================
  # 3b. PCA with Significant DMLs only (Top hits)
  # ==============================================================================
  logit("  3b. PCA with Significant DMLs...")
  
  sig_dmls_pca <- dml_res[pval < DSS_P_VAL & abs(diff) > DSS_DELTA]
  sig_dmls_pca <- sig_dmls_pca[order(pval)]
  
  # Top 100 또는 전체 유의한 DML 중 작은 수
  top_n_pca <- min(100, nrow(sig_dmls_pca))
  
  if(top_n_pca > 10) {
    top_dmls_pca <- sig_dmls_pca[1:top_n_pca]
    
    gr_bs <- granges(BSobj)
    top_gr_pca <- GRanges(seqnames = top_dmls_pca$chr, 
                          ranges = IRanges(start = top_dmls_pca$pos, width = 1))
    
    overlaps_pca <- findOverlaps(top_gr_pca, gr_bs)
    matched_idx_pca <- subjectHits(overlaps_pca)
    
    if(length(matched_idx_pca) > 10) {
      M_mat_sig <- as.matrix(getMeth(BSobj[matched_idx_pca, ], type = "raw"))
      M_mat_sig[is.na(M_mat_sig)] <- 0.5
      
      row_vars_sig <- apply(M_mat_sig, 1, var, na.rm = TRUE)
      M_mat_sig <- M_mat_sig[row_vars_sig > 0 & !is.na(row_vars_sig), , drop = FALSE]
      
      if(nrow(M_mat_sig) > 10) {
        pca_res_sig <- prcomp(t(M_mat_sig), scale. = TRUE)
        pca_df_sig <- data.frame(
          PC1 = pca_res_sig$x[, 1],
          PC2 = pca_res_sig$x[, 2],
          Sample = colnames(M_mat_sig),
          Group = targets$Group
        )
        
        var_explained_sig <- round(100 * pca_res_sig$sdev^2 / sum(pca_res_sig$sdev^2), 1)
        
        p_pca_sig <- ggplot(pca_df_sig, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
          geom_point(size = 4) +
          geom_text(vjust = -1, size = 3) +
          scale_color_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
          labs(title = paste0("PCA: Top ", nrow(M_mat_sig), " Significant DMLs Only"),
               subtitle = paste0("(p < ", DSS_P_VAL, ", |diff| > ", DSS_DELTA, ")"),
               x = paste0("PC1 (", var_explained_sig[1], "%)"),
               y = paste0("PC2 (", var_explained_sig[2], "%)")) +
          theme_bw() +
          theme(legend.position = "bottom")
        
        ggsave(file.path(plot_dir, "03b_PCA_SignificantDMLs.png"), p_pca_sig, width = 10, height = 8, dpi = 300)
        logit("     -> PCA with ", nrow(M_mat_sig), " significant DMLs created")
      }
    }
  }
  
  # ==============================================================================
  # 4. Top DML Heatmap
  # ==============================================================================
  logit("  4. Top DML Heatmap...")
  
  sig_dmls <- dml_res[pval < DSS_P_VAL & abs(diff) > DSS_DELTA]
  sig_dmls <- sig_dmls[order(pval)]
  
  top_n_dml <- min(100, nrow(sig_dmls))
  
  if(top_n_dml > 5) {
    top_dmls <- sig_dmls[1:top_n_dml]
    
    gr_bs <- granges(BSobj)
    top_gr <- GRanges(seqnames = top_dmls$chr, 
                      ranges = IRanges(start = top_dmls$pos, width = 1))
    
    overlaps <- findOverlaps(top_gr, gr_bs)
    matched_idx <- subjectHits(overlaps)
    
    if(length(matched_idx) > 5) {
      dml_meth <- as.matrix(getMeth(BSobj[matched_idx, ], type = "raw"))
      rownames(dml_meth) <- paste0(top_dmls$chr[queryHits(overlaps)], ":", 
                                    top_dmls$pos[queryHits(overlaps)])
      
      dml_meth[is.na(dml_meth)] <- 0.5
      
      row_vars <- apply(dml_meth, 1, var, na.rm = TRUE)
      dml_meth <- dml_meth[row_vars > 0 & !is.na(row_vars), , drop = FALSE]
      
      if(nrow(dml_meth) > 5) {
        sample_order <- order(targets$Group, decreasing = FALSE)
        dml_meth_ordered <- dml_meth[, sample_order, drop = FALSE]
        targets_ordered <- targets[sample_order, ]
        
        anno_col <- data.frame(Group = targets_ordered$Group, row.names = targets_ordered$SampleID)
        anno_colors <- list(Group = c("Control" = "#4DAF4A", "Case" = "#E41A1C"))
        
        n_case <- sum(targets_ordered$Group == "Case")
        
        smooth_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                            "#F7F7F7", 
                                            "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(256)
        
        png(file.path(plot_dir, "04_Heatmap_TopDMLs.png"), width = 1400, height = 1000, res = 150)
        pheatmap(dml_meth_ordered,
                 scale = "row",
                 clustering_distance_rows = "correlation",
                 cluster_cols = FALSE,
                 annotation_col = anno_col,
                 annotation_colors = anno_colors,
                 color = smooth_colors,
                 main = paste("Top", nrow(dml_meth), "Significant DMLs (Case vs Control)"),
                 fontsize_row = 5,
                 fontsize_col = 7,
                 show_rownames = FALSE,
                 gaps_col = n_case,
                 border_color = NA)
        dev.off()
      }
    }
  }
  
  # ==============================================================================
  # 5. DMR Heatmaps (04b and 04c)
  # ==============================================================================
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  5. DMR Heatmaps...")
    
    dmr_gr <- GRanges(seqnames = dmrs$chr, 
                      ranges = IRanges(start = dmrs$start, end = dmrs$end))
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
    gene_names <- character(nrow(dmrs))
    for(i in 1:nrow(dmrs)) {
      gene_name <- get_gene_at_position(dmrs$chr[i], dmrs$start[i], dmrs$end[i])
      gene_names[i] <- if(!is.na(gene_name) && gene_name != "") gene_name else paste0(dmrs$chr[i], ":", dmrs$start[i])
      Sys.sleep(0.1)
    }
    
    if(any(duplicated(gene_names))) {
      dup_idx <- which(duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE))
      for(idx in dup_idx) {
        gene_names[idx] <- paste0(gene_names[idx], " (", dmrs$chr[idx], ":", dmrs$start[idx], ")")
      }
    }
    
    rownames(dmr_meth_mat) <- gene_names
    colnames(dmr_meth_mat) <- colnames(BSobj)
    
    for(i in 1:nrow(dmr_meth_mat)) {
      na_idx <- is.na(dmr_meth_mat[i, ])
      if(any(na_idx)) {
        row_mean <- mean(dmr_meth_mat[i, !na_idx], na.rm = TRUE)
        if(is.na(row_mean)) row_mean <- 0.5
        dmr_meth_mat[i, na_idx] <- row_mean
      }
    }
    dmr_meth_mat[is.na(dmr_meth_mat)] <- 0.5
    
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
      
      row_vars <- apply(dmr_meth_ordered, 1, var, na.rm = TRUE)
      do_cluster_rows <- all(row_vars > 1e-10, na.rm = TRUE) && nrow(dmr_meth_ordered) > 2
      
      # 04b: Standard DMR heatmap (with clustering)
      png(file.path(plot_dir, "04b_Heatmap_DMR_Methylation.png"), 
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
               main = paste("DMR Methylation Levels (n=", nrow(dmr_meth_ordered), " DMRs)"),
               fontsize_row = 8,
               fontsize_col = 6,
               show_rownames = TRUE,
               show_colnames = FALSE,
               gaps_col = n_case,
               border_color = NA,
               legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
               legend_labels = c("0", "0.25", "0.5", "0.75", "1"))
      dev.off()
      
      # 04c: Grouped DMR heatmap (Hyper/Hypo separated)
      # DSS에서 diff.Methy = meanMethy1(Control) - meanMethy2(Case)
      # diff.Methy < 0: Case > Control → Hypermethylation
      dmr_direction <- ifelse(dmrs$diff.Methy < 0, "Hypermethylation", "Hypomethylation")
      
      hyper_idx <- which(dmr_direction == "Hypermethylation")
      hypo_idx <- which(dmr_direction == "Hypomethylation")
      hyper_order <- hyper_idx[order(-abs(dmrs$diff.Methy[hyper_idx]))]
      hypo_order <- hypo_idx[order(-abs(dmrs$diff.Methy[hypo_idx]))]
      row_order <- c(hyper_order, hypo_order)
      
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
      
      png(file.path(plot_dir, "04c_Heatmap_DMR_Grouped.png"), 
          width = 1600, height = heatmap_height, res = 150)
      pheatmap(dmr_meth_grouped,
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_col = anno_col,
               annotation_row = anno_row,
               annotation_colors = anno_colors_grouped,
               color = viridis_colors,
               main = paste0("DMR Methylation Levels - Grouped by Direction\n",
                            "(Hypermethylation: n=", n_hyper, 
                            ", Hypomethylation: n=", length(hypo_order), ")"),
               fontsize_row = 8,
               fontsize_col = 6,
               show_rownames = TRUE,
               show_colnames = FALSE,
               gaps_col = n_case,
               gaps_row = n_hyper,
               border_color = NA,
               legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
               legend_labels = c("0", "0.25", "0.5", "0.75", "1"))
      dev.off()
      
      logit("     -> Heatmaps created: ", n_hyper, " Hyper + ", length(hypo_order), " Hypo DMRs")
      
      # ==============================================================
      # 5b. PCA with DMR methylation values only
      # ==============================================================
      logit("  5b. PCA with DMR regions only...")
      
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
        labs(title = paste0("PCA: ", nrow(dmr_meth_ordered), " DMR Regions Only"),
             subtitle = paste0("(Hypermethylation: ", n_hyper, ", Hypomethylation: ", length(hypo_order), ")"),
             x = paste0("PC1 (", var_explained_dmr[1], "%)"),
             y = paste0("PC2 (", var_explained_dmr[2], "%)")) +
        theme_bw() +
        theme(legend.position = "bottom")
      
      ggsave(file.path(plot_dir, "03c_PCA_DMRs.png"), p_pca_dmr, width = 10, height = 8, dpi = 300)
      logit("     -> PCA with ", nrow(dmr_meth_ordered), " DMRs created")
    }
  }
  
  # ==============================================================================
  # 6. DMR Length Distribution
  # ==============================================================================
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  6. DMR Length Distribution...")
    
    dmrs[, length := end - start]
    
    p_length <- ggplot(dmrs, aes(x = length)) +
      geom_histogram(bins = 50, fill = "#3288BD", color = "white", alpha = 0.8) +
      scale_x_log10() +
      labs(title = "DMR Length Distribution",
           x = "DMR Length (bp, log scale)",
           y = "Count") +
      theme_bw()
    
    ggsave(file.path(plot_dir, "05_DMR_Length_Distribution.png"), p_length, width = 8, height = 6, dpi = 300)
  }
  
  # ==============================================================================
  # 7. Chromosome Distribution
  # ==============================================================================
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  7. Chromosome Distribution...")
    
    chr_counts <- dmrs[, .N, by = chr]
    chr_counts[, chr_num := gsub("chr", "", chr)]
    chr_counts[, chr_num := factor(chr_num, levels = c(1:22, "X", "Y"))]
    chr_counts <- chr_counts[!is.na(chr_num)]
    
    p_chr <- ggplot(chr_counts, aes(x = chr_num, y = N, fill = N)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#FEE08B", high = "#D53E4F") +
      labs(title = "DMR Distribution by Chromosome",
           x = "Chromosome",
           y = "Number of DMRs") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(plot_dir, "06_DMR_Chromosome_Distribution.png"), p_chr, width = 10, height = 6, dpi = 300)
  }
  
  # ==============================================================================
  # 8. Hyper/Hypo Pie Chart
  # ==============================================================================
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  8. Hyper/Hypo Pie Chart...")
    
    # DSS에서 diff.Methy = meanMethy1(Control) - meanMethy2(Case)
    dmrs[, direction := ifelse(diff.Methy < 0, "Hypermethylated", "Hypomethylated")]
    direction_counts <- dmrs[, .N, by = direction]
    direction_counts[, pct := round(N / sum(N) * 100, 1)]
    direction_counts[, label := paste0(direction, "\n(n=", N, ", ", pct, "%)")]
    
    p_pie <- ggplot(direction_counts, aes(x = "", y = N, fill = direction)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = c("Hypermethylated" = "#E41A1C", "Hypomethylated" = "#377EB8")) +
      geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
      labs(title = "DMR Direction: Hyper vs Hypo Methylation (in Case)") +
      theme_void() +
      theme(legend.position = "none")
    
    ggsave(file.path(plot_dir, "07_DMR_Direction_PieChart.png"), p_pie, width = 8, height = 8, dpi = 300)
  }
  
  # ==============================================================================
  # 9. Top DMR Detail Plots (PDF)
  # ==============================================================================
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  9. Top DMR Detail Plots...")
    
    dml_df <- as.data.frame(dml_res)
    dmrs_df <- as.data.frame(dmrs)
    dmrs_df <- dmrs_df[order(-abs(dmrs_df$diff.Methy)), ]
    
    top_n <- min(10, nrow(dmrs_df))
    
    pdf(file.path(plot_dir, "08_Top_DMRs_Detail.pdf"), width = 12, height = 6)
    
    for(i in 1:top_n) {
      dmr <- dmrs_df[i, ]
      
      tryCatch({
        showOneDMR(dmrs_df, dml_df, i,
                   ext = 500,
                   main = paste0("Rank ", i, ": ", dmr$chr, ":", dmr$start, "-", dmr$end,
                                " (diff=", round(dmr$diff.Methy, 3), ")"))
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error plotting DMR", i))
      })
    }
    
    dev.off()
  }
  
  # ==============================================================================
  # 10. Summary Statistics
  # ==============================================================================
  logit("  10. Summary Statistics...")
  
  # DSS에서 diff = meanMethy1(Control) - meanMethy2(Case)
  summary_stats <- data.frame(
    Metric = c(
      "Total CpG Sites Tested",
      "Significant DMLs (p < threshold)",
      "Hypermethylated DMLs (Case > Control)",
      "Hypomethylated DMLs (Case < Control)",
      "Total DMRs",
      "Hypermethylated DMRs (Case > Control)",
      "Hypomethylated DMRs (Case < Control)",
      "Median DMR Length (bp)",
      "Control Samples",
      "Case Samples"
    ),
    Value = c(
      nrow(dml_res),
      nrow(dml_res[pval < DSS_P_VAL & abs(diff) > DSS_DELTA]),
      nrow(dml_res[pval < DSS_P_VAL & diff < -DSS_DELTA]),
      nrow(dml_res[pval < DSS_P_VAL & diff > DSS_DELTA]),
      if(!is.null(dmrs)) nrow(dmrs) else 0,
      if(!is.null(dmrs)) nrow(dmrs[diff.Methy < 0]) else 0,
      if(!is.null(dmrs)) nrow(dmrs[diff.Methy > 0]) else 0,
      if(!is.null(dmrs) && nrow(dmrs) > 0) median(dmrs$end - dmrs$start) else NA,
      sum(targets$Group == "Control"),
      sum(targets$Group == "Case")
    )
  )
  
  fwrite(summary_stats, file.path(plot_dir, "00_Summary_Statistics.csv"))
  
  logit("  -> Done!")
}

logit("========================================")
logit("All folders processed successfully!")
