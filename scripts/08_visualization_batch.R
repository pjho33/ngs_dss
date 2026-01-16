# ==============================================================================
# Script: 08_visualization_batch.R
# Purpose: Run FULL visualization for all result subfolders
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
  library(GenomicRanges)
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
subdirs <- subdirs[grepl("^p[0-9]", basename(subdirs))]  # p로 시작하는 폴더만

logit("Found", length(subdirs), "result folders to process")

for(subdir in subdirs) {
  logit("========================================")
  logit("Processing:", basename(subdir))
  
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
  
  logit("  Loaded:", nrow(dml_res), "DMLs,", if(!is.null(dmrs)) nrow(dmrs) else 0, "DMRs")
  
  # Parameters from folder name
  folder_params <- strsplit(basename(subdir), "_")[[1]]
  DSS_P_VAL <- as.numeric(gsub("p", "", folder_params[1]))
  DSS_DELTA <- as.numeric(gsub("d", "", folder_params[2]))
  
  # ------------------------------------------------------------------------------
  # 1. Volcano Plot
  # ------------------------------------------------------------------------------
  logit("  1. Volcano Plot...")
  
  # DSS에서 diff = meanMethy1(Control) - meanMethy2(Case)
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
  
  # ------------------------------------------------------------------------------
  # 2. Hyper/Hypo Pie Chart (DMR)
  # ------------------------------------------------------------------------------
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  2. Pie Chart...")
    
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
  
  # ------------------------------------------------------------------------------
  # 3. DMR Grouped Heatmap
  # ------------------------------------------------------------------------------
  if(!is.null(dmrs) && nrow(dmrs) > 0) {
    logit("  3. Grouped DMR Heatmap...")
    
    # DMR 영역에 대한 GRanges 생성
    dmr_gr <- GRanges(seqnames = dmrs$chr, 
                      ranges = IRanges(start = dmrs$start, end = dmrs$end))
    bs_gr <- granges(BSobj)
    
    # 각 DMR 영역 내 CpG들의 평균 메틸레이션 계산
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
    
    # Gene annotation (Ensembl REST API)
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
    
    gene_names <- character(nrow(dmrs))
    for(i in 1:nrow(dmrs)) {
      gene_name <- get_gene_at_position(dmrs$chr[i], dmrs$start[i], dmrs$end[i])
      gene_names[i] <- if(!is.na(gene_name) && gene_name != "") gene_name else paste0(dmrs$chr[i], ":", dmrs$start[i])
      Sys.sleep(0.1)
    }
    
    # 중복 이름 처리
    if(any(duplicated(gene_names))) {
      dup_idx <- which(duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE))
      for(idx in dup_idx) {
        gene_names[idx] <- paste0(gene_names[idx], " (", dmrs$chr[idx], ":", dmrs$start[idx], ")")
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
    
    if(nrow(dmr_meth_mat) > 1) {
      # 샘플 정렬
      sample_order <- order(targets$Group, decreasing = FALSE)
      dmr_meth_ordered <- dmr_meth_mat[, sample_order, drop = FALSE]
      targets_ordered <- targets[sample_order, ]
      
      # DMR 방향 정보 (수정된 버전)
      # diff.Methy < 0: Case > Control → Hypermethylation
      dmr_direction <- ifelse(dmrs$diff.Methy < 0, "Hypermethylation", "Hypomethylation")
      
      hyper_idx <- which(dmr_direction == "Hypermethylation")
      hypo_idx <- which(dmr_direction == "Hypomethylation")
      hyper_order <- hyper_idx[order(-abs(dmrs$diff.Methy[hyper_idx]))]
      hypo_order <- hypo_idx[order(-abs(dmrs$diff.Methy[hypo_idx]))]
      row_order <- c(hyper_order, hypo_order)
      
      dmr_meth_grouped <- dmr_meth_ordered[row_order, , drop = FALSE]
      dmr_direction_ordered <- dmr_direction[row_order]
      
      anno_col <- data.frame(Group = targets_ordered$Group, row.names = targets_ordered$SampleID)
      anno_row <- data.frame(
        Direction = factor(dmr_direction_ordered, levels = c("Hypermethylation", "Hypomethylation")),
        row.names = rownames(dmr_meth_grouped)
      )
      anno_colors <- list(
        Group = c("Control" = "#4DAF4A", "Case" = "#E41A1C"),
        Direction = c("Hypermethylation" = "#D73027", "Hypomethylation" = "#4575B4")
      )
      
      n_case <- sum(targets_ordered$Group == "Case")
      n_hyper <- length(hyper_order)
      
      viridis_colors <- colorRampPalette(c("#440154", "#482878", "#3E4A89", "#31688E", 
                                            "#26828E", "#1F9E89", "#35B779", "#6DCD59", 
                                            "#B4DE2C", "#FDE725"))(256)
      
      heatmap_height <- max(600, min(1500, 100 + nrow(dmr_meth_grouped) * 25))
      
      png(file.path(plot_dir, "04c_Heatmap_DMR_Grouped.png"), 
          width = 1600, height = heatmap_height, res = 150)
      pheatmap(dmr_meth_grouped,
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_col = anno_col,
               annotation_row = anno_row,
               annotation_colors = anno_colors,
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
      logit("     -> Grouped Heatmap:", n_hyper, "Hyper +", length(hypo_order), "Hypo")
    }
  }
  
  # ------------------------------------------------------------------------------
  # 4. Summary Statistics
  # ------------------------------------------------------------------------------
  logit("  4. Summary Statistics...")
  
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
