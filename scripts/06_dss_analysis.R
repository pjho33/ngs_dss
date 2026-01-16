# ==============================================================================
# Script: 06_dss_analysis.R
# Purpose: DSS analysis using patients/normals folder structure
# Input: results/methylation/patients/ and results/methylation/normals/
# Output: results/dss_results/
# ==============================================================================

suppressPackageStartupMessages({
  library(DSS)
  library(bsseq)
  library(data.table)
  library(GenomicRanges)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# 경로 설정
BASE_DIR <- "/home/pjho3/projects/ngs_dss"
PATIENTS_DIR <- file.path(BASE_DIR, "results/methylation/patients")
NORMALS_DIR <- file.path(BASE_DIR, "results/methylation/normals")
OUT_DIR <- file.path(BASE_DIR, "results/dss_results")

# DSS 파라미터
DSS_SMOOTHING <- TRUE
DSS_DELTA <- 0.1
DSS_P_VAL <- 0.001
DMR_MIN_LEN <- 50
DMR_MIN_CG <- 3
DMR_DIS_MERGE <- 100

# 표준 염색체만 사용
STD_CHRS <- c(paste0("chr", c(1:22, "X", "Y")), as.character(c(1:22, "X", "Y")))

# ------------------------------------------------------------------------------
# 2. Load Coverage Files
# ------------------------------------------------------------------------------
logit("=== DSS Analysis Pipeline ===")
logit("Loading coverage files from patients/normals folders...")

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 파일 목록 가져오기
patient_files <- list.files(PATIENTS_DIR, pattern = "\\.cov\\.gz$", full.names = TRUE)
normal_files <- list.files(NORMALS_DIR, pattern = "\\.cov\\.gz$", full.names = TRUE)

logit(paste("Patient files:", length(patient_files)))
logit(paste("Normal files:", length(normal_files)))

if(length(patient_files) == 0 || length(normal_files) == 0) {
  stop("No coverage files found in patients or normals folder!")
}

# ------------------------------------------------------------------------------
# 3. Read and Process Data
# ------------------------------------------------------------------------------
read_bismark_cov <- function(file_path, sample_name) {
  logit(paste("  Reading:", sample_name))
  dt <- fread(file_path, header = FALSE, 
              col.names = c("chr", "start", "end", "meth_pct", "count_m", "count_u"))
  dt[, chr := as.character(chr)]
  dt <- dt[chr %in% STD_CHRS]
  dt[, N := count_m + count_u]
  dt[, X := count_m]
  return(dt[, .(chr, pos = start, N, X)])
}

# 샘플 이름 추출 함수
extract_sample_name <- function(filepath) {
  basename(filepath) |> 
    gsub("_L005.*", "", x = _) |>
    gsub("_S[0-9]+.*", "", x = _)
}

logit("Reading patient samples...")
patient_data <- list()
patient_names <- c()
for(f in patient_files) {
  sname <- extract_sample_name(f)
  patient_names <- c(patient_names, sname)
  patient_data[[sname]] <- read_bismark_cov(f, sname)
}

logit("Reading normal samples...")
normal_data <- list()
normal_names <- c()
for(f in normal_files) {
  sname <- extract_sample_name(f)
  normal_names <- c(normal_names, sname)
  normal_data[[sname]] <- read_bismark_cov(f, sname)
}

# ------------------------------------------------------------------------------
# 4. Create BSseq Object
# ------------------------------------------------------------------------------
logit("Creating BSseq objects...")

# DSS용 BSseq 객체 생성 함수
create_bsseq_from_list <- function(data_list, sample_names) {
  all_data <- rbindlist(data_list, idcol = "sample")
  
  # 모든 위치 추출
  all_pos <- unique(all_data[, .(chr, pos)])
  setkey(all_pos, chr, pos)
  
  # M과 Cov 매트릭스 생성
  n_sites <- nrow(all_pos)
  n_samples <- length(sample_names)
  
  M_mat <- matrix(0L, nrow = n_sites, ncol = n_samples)
  Cov_mat <- matrix(0L, nrow = n_sites, ncol = n_samples)
  colnames(M_mat) <- sample_names
  colnames(Cov_mat) <- sample_names
  
  for(i in seq_along(sample_names)) {
    sname <- sample_names[i]
    dt <- data_list[[sname]]
    setkey(dt, chr, pos)
    
    # 매칭
    idx <- all_pos[dt, on = .(chr, pos), which = TRUE, nomatch = 0]
    matched <- dt[all_pos, on = .(chr, pos), nomatch = NA]
    
    M_mat[, i] <- ifelse(is.na(matched$X), 0L, as.integer(matched$X))
    Cov_mat[, i] <- ifelse(is.na(matched$N), 0L, as.integer(matched$N))
  }
  
  # GRanges 생성
  gr <- GRanges(seqnames = all_pos$chr, 
                ranges = IRanges(start = all_pos$pos, width = 1))
  
  BSseq(M = M_mat, Cov = Cov_mat, gr = gr, sampleNames = sample_names)
}

# 전체 데이터 합치기
all_data <- c(patient_data, normal_data)
all_names <- c(patient_names, normal_names)

logit("Building combined BSseq object (this may take a while)...")
bs_obj <- create_bsseq_from_list(all_data, all_names)

# 그룹 정보 추가
pData(bs_obj)$Group <- c(rep("Case", length(patient_names)), 
                          rep("Control", length(normal_names)))

logit(paste("BSseq object created:", length(bs_obj), "CpGs,", ncol(bs_obj), "samples"))

# 메모리 정리
rm(all_data, patient_data, normal_data)
gc()

# ------------------------------------------------------------------------------
# 5. Run DSS Analysis
# ------------------------------------------------------------------------------
logit("=== Running DSS Analysis ===")
logit(paste("Parameters: smoothing=", DSS_SMOOTHING, ", delta=", DSS_DELTA, ", p.threshold=", DSS_P_VAL))

logit("1. Running DML Test (Case vs Control)...")
logit(paste("   Case samples:", length(patient_names)))
logit(paste("   Control samples:", length(normal_names)))

dml_res <- DMLtest(bs_obj, group1 = normal_names, group2 = patient_names, smoothing = DSS_SMOOTHING)

# 전체 DML 결과 저장
fwrite(as.data.table(dml_res), file.path(OUT_DIR, "All_DMLs.tsv.gz"), sep = "\t", compress = "gzip")
logit(paste("   -> Total Sites tested:", nrow(dml_res)))

logit("2. Filtering Significant DMLs...")
sig_dmls <- dml_res[!is.na(dml_res$pval) & dml_res$pval < DSS_P_VAL & abs(dml_res$diff) > DSS_DELTA, ]
fwrite(sig_dmls, file.path(OUT_DIR, "Significant_DMLs.tsv"), sep = "\t")
logit(paste("   -> Significant DMLs:", nrow(sig_dmls)))

logit("3. Calling DMRs...")
dmrs <- callDMR(dml_res, p.threshold = DSS_P_VAL, delta = DSS_DELTA,
                minlen = DMR_MIN_LEN, minCG = DMR_MIN_CG, dis.merge = DMR_DIS_MERGE)

if(!is.null(dmrs) && is.data.frame(dmrs) && nrow(dmrs) > 0) {
  fwrite(dmrs, file.path(OUT_DIR, "Final_DMRs.tsv"), sep = "\t")
  logit(paste("   -> DMRs found:", nrow(dmrs)))
  
  # BED 파일 생성 (시각화용)
  bed_df <- as.data.frame(dmrs[, c("chr", "start", "end", "diff.Methy")])
  bed_df$name <- paste0("DMR_", 1:nrow(bed_df))
  bed_df$score <- round(bed_df$diff.Methy * 1000)
  write.table(bed_df[, c(1, 2, 3, 5, 6)], file.path(OUT_DIR, "DMRs_visualization.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  logit("   -> BED file created for visualization")
} else {
  logit("   -> No significant DMRs found with current parameters")
}

# 샘플 정보 저장
sample_info <- data.table(
  SampleID = all_names,
  Group = c(rep("Case", length(patient_names)), rep("Control", length(normal_names))),
  Source = c(rep("patients", length(patient_names)), rep("normals", length(normal_names)))
)
fwrite(sample_info, file.path(OUT_DIR, "sample_info.tsv"), sep = "\t")

logit("======================================================")
logit(" SUCCESS! DSS Analysis Complete")
logit(paste(" Results saved to:", OUT_DIR))
logit("======================================================")
