# ==============================================================================
# Script: 06_import_bsseq.R (Optimized Version)
# Purpose: Create BSseq object using read.bismark for memory-efficient DSS analysis
# Output: 'r_objects/dss_h5_store/bsseq.rds'
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(bsseq)
  library(data.table)
  library(GenomicRanges)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

# 입력/출력 경로 설정 (patients/normals 폴더 구조 사용)
patients_dir <- file.path(paths$output_dir, "methylation", "patients")
normals_dir <- file.path(paths$output_dir, "methylation", "normals")
OUT_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# ------------------------------------------------------------------------------
# 2. Setup & Validation
# ------------------------------------------------------------------------------
logit("=== [Step 1] Initializing BSseq Data Pipeline ===")

if(!dir.exists(patients_dir)) stop(paste("Directory not found:", patients_dir))
if(!dir.exists(normals_dir)) stop(paste("Directory not found:", normals_dir))

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 파일 리스트 확보 (patients/normals 폴더에서)
patient_files <- list.files(patients_dir, full.names = TRUE, 
                            pattern = "\\.cov$|\\.cov\\.gz$")
normal_files <- list.files(normals_dir, full.names = TRUE, 
                           pattern = "\\.cov$|\\.cov\\.gz$")

if(length(patient_files) == 0) stop(paste("No coverage files found in:", patients_dir))
if(length(normal_files) == 0) stop(paste("No coverage files found in:", normals_dir))

# 샘플 이름 추출
patient_names <- gsub("_L005.*", "", gsub("_S[0-9]+.*", "", basename(patient_files)))
normal_names <- gsub("_L005.*", "", gsub("_S[0-9]+.*", "", basename(normal_files)))

logit(paste("Detected Samples:", length(patient_files) + length(normal_files)))
logit(paste("  - Patients (Case):", length(patient_files)))
logit(paste("  - Normals (Control):", length(normal_files)))

# ------------------------------------------------------------------------------
# 3. Load Data via Optimized bsseq::read.bismark
# ------------------------------------------------------------------------------
logit("=== [Step 2] Loading Data via read.bismark ===")

all_files <- c(patient_files, normal_files)
all_names <- c(patient_names, normal_names)

# read.bismark 함수 사용 (백엔드 최적화로 속도가 훨씬 빠름)
bs_obj <- read.bismark(
  files = all_files,
  rmZeroCov = TRUE,      # 커버리지 0인 지역 제거 (속도 향상)
  strandCollapse = TRUE, # strand 합치기
  verbose = TRUE
)

# 샘플 이름 설정
sampleNames(bs_obj) <- all_names

# 그룹 정보 추가
pData(bs_obj)$Group <- c(rep("Case", length(patient_names)), 
                         rep("Control", length(normal_names)))

# ------------------------------------------------------------------------------
# 4. Filter Standard Chromosomes
# ------------------------------------------------------------------------------
logit("=== [Step 3] Filtering Standard Chromosomes ===")

std_chrs <- c(paste0("chr", c(1:22, "X", "Y")), as.character(c(1:22, "X", "Y")))
bs_obj <- bs_obj[seqnames(bs_obj) %in% std_chrs, ]

logit(paste("BSseq object created:", nrow(bs_obj), "CpGs,", ncol(bs_obj), "samples"))

# ------------------------------------------------------------------------------
# 5. Save Results
# ------------------------------------------------------------------------------
logit("=== [Step 4] Saving BSseq Object ===")

saveRDS(bs_obj, file.path(OUT_DIR, "bsseq.rds"))

# targets 정보도 저장 (07번 스크립트에서 사용)
targets <- data.table(
  SampleID = all_names,
  Group = c(rep("Case", length(patient_names)), rep("Control", length(normal_names))),
  FilePath = all_files
)
saveRDS(targets, file.path(OUT_DIR, "targets.rds"))

logit("======================================================")
logit(" SUCCESS! BSseq object created at: ", OUT_DIR)
logit(paste(" Total CpGs:", format(nrow(bs_obj), big.mark=",")))
logit(paste(" Total Samples:", ncol(bs_obj)))
logit(" You can now run 07_dss_modeling.R for analysis.")
logit("======================================================")
