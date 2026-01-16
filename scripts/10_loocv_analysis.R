# ==============================================================================
# Script: 10_loocv_analysis.R
# Purpose: LOOCV (Leave-One-Out Cross-Validation) for model stability assessment
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(bsseq)
  library(data.table)
  library(pROC)
  library(GenomicRanges)
})

paths <- read_yaml("config/paths.yaml")
H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
DMRSEQ_DIR <- file.path(paths$output_dir, "dss_results", "dmrseq", "20251228_cut0.05_minR3_bp1000")

# Load data
BSobj <- readRDS(file.path(H5_DIR, "bsseq.rds"))
targets <- readRDS(file.path(H5_DIR, "targets.rds"))
sens_df <- fread(file.path(DMRSEQ_DIR, "sensitivity_analysis.tsv"))

bs_gr <- granges(BSobj)

get_meth_data <- function(chr, start, end) {
  dmr_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  overlaps <- findOverlaps(dmr_gr, bs_gr)
  cpg_idx <- subjectHits(overlaps)
  if(length(cpg_idx) > 0) {
    meth_vals <- as.matrix(getMeth(BSobj[cpg_idx, ], type = "raw"))
    return(colMeans(meth_vals, na.rm = TRUE))
  }
  return(NULL)
}

# Get BNC1 and NXPH1 data
bnc1_row <- sens_df[Gene == "BNC1"]
nxph1_row <- sens_df[Gene == "NXPH1"]

bnc1_meth <- get_meth_data(bnc1_row$Chr, bnc1_row$Start, bnc1_row$End)
nxph1_meth <- get_meth_data(nxph1_row$Chr, nxph1_row$Start, nxph1_row$End)

# Create data frame
model_data <- data.frame(
  BNC1 = bnc1_meth,
  NXPH1 = nxph1_meth,
  Group = ifelse(targets$Group == "Case", 1, 0)
)
model_data <- model_data[complete.cases(model_data), ]

n <- nrow(model_data)
cat("=== LOOCV (Leave-One-Out Cross-Validation) ===\n\n")
cat(sprintf("Total samples: %d (Case: %d, Control: %d)\n\n", n, sum(model_data$Group), sum(model_data$Group == 0)))

# LOOCV for each model
loocv_results <- data.frame(
  actual = model_data$Group,
  pred_bnc1 = numeric(n),
  pred_nxph1 = numeric(n),
  pred_combined = numeric(n)
)

cat("Running LOOCV (", n, " iterations)...\n", sep = "")

for(i in 1:n) {
  # Leave one out
  train <- model_data[-i, ]
  test <- model_data[i, ]
  
  # Fit models on training data
  suppressWarnings({
    model_bnc1 <- glm(Group ~ BNC1, data = train, family = binomial)
    model_nxph1 <- glm(Group ~ NXPH1, data = train, family = binomial)
    model_combined <- glm(Group ~ BNC1 + NXPH1, data = train, family = binomial)
  })
  
  # Predict on left-out sample
  loocv_results$pred_bnc1[i] <- predict(model_bnc1, newdata = test, type = "response")
  loocv_results$pred_nxph1[i] <- predict(model_nxph1, newdata = test, type = "response")
  loocv_results$pred_combined[i] <- predict(model_combined, newdata = test, type = "response")
}

cat("Done!\n\n")

# Calculate LOOCV AUC
roc_bnc1_cv <- roc(loocv_results$actual, loocv_results$pred_bnc1, quiet = TRUE)
roc_nxph1_cv <- roc(loocv_results$actual, loocv_results$pred_nxph1, quiet = TRUE)
roc_combined_cv <- roc(loocv_results$actual, loocv_results$pred_combined, quiet = TRUE)

cat("=== LOOCV AUC Results ===\n")
cat(sprintf("BNC1 only:      AUC = %.4f\n", auc(roc_bnc1_cv)))
cat(sprintf("NXPH1 only:     AUC = %.4f\n", auc(roc_nxph1_cv)))
cat(sprintf("Combined:       AUC = %.4f\n", auc(roc_combined_cv)))

# Calculate accuracy, sensitivity, specificity at optimal threshold
get_metrics <- function(actual, predicted) {
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  
  pred_class <- ifelse(predicted >= coords$threshold, 1, 0)
  accuracy <- mean(pred_class == actual)
  
  tp <- sum(pred_class == 1 & actual == 1)
  tn <- sum(pred_class == 0 & actual == 0)
  fp <- sum(pred_class == 1 & actual == 0)
  fn <- sum(pred_class == 0 & actual == 1)
  
  list(
    AUC = auc(roc_obj),
    Accuracy = accuracy,
    Sensitivity = coords$sensitivity,
    Specificity = coords$specificity,
    PPV = tp / (tp + fp),
    NPV = tn / (tn + fn)
  )
}

metrics_bnc1 <- get_metrics(loocv_results$actual, loocv_results$pred_bnc1)
metrics_nxph1 <- get_metrics(loocv_results$actual, loocv_results$pred_nxph1)
metrics_combined <- get_metrics(loocv_results$actual, loocv_results$pred_combined)

cat("\n=== LOOCV Performance Metrics ===\n\n")
cat(sprintf("%-20s %12s %12s %12s\n", "Metric", "BNC1", "NXPH1", "Combined"))
cat(paste(rep("-", 60), collapse = ""), "\n")
cat(sprintf("%-20s %11.1f%% %11.1f%% %11.1f%%\n", "Accuracy", metrics_bnc1$Accuracy*100, metrics_nxph1$Accuracy*100, metrics_combined$Accuracy*100))
cat(sprintf("%-20s %11.1f%% %11.1f%% %11.1f%%\n", "Sensitivity", metrics_bnc1$Sensitivity*100, metrics_nxph1$Sensitivity*100, metrics_combined$Sensitivity*100))
cat(sprintf("%-20s %11.1f%% %11.1f%% %11.1f%%\n", "Specificity", metrics_bnc1$Specificity*100, metrics_nxph1$Specificity*100, metrics_combined$Specificity*100))
cat(sprintf("%-20s %11.1f%% %11.1f%% %11.1f%%\n", "PPV", metrics_bnc1$PPV*100, metrics_nxph1$PPV*100, metrics_combined$PPV*100))
cat(sprintf("%-20s %11.1f%% %11.1f%% %11.1f%%\n", "NPV", metrics_bnc1$NPV*100, metrics_nxph1$NPV*100, metrics_combined$NPV*100))
cat(sprintf("%-20s %12.4f %12.4f %12.4f\n", "AUC", metrics_bnc1$AUC, metrics_nxph1$AUC, metrics_combined$AUC))

# Compare with original (non-CV) AUC
cat("\n=== Model Stability (Original vs LOOCV) ===\n")

suppressWarnings({
  full_model_bnc1 <- glm(Group ~ BNC1, data = model_data, family = binomial)
  full_model_nxph1 <- glm(Group ~ NXPH1, data = model_data, family = binomial)
  full_model_combined <- glm(Group ~ BNC1 + NXPH1, data = model_data, family = binomial)
})

orig_auc_bnc1 <- auc(roc(model_data$Group, predict(full_model_bnc1, type="response"), quiet=TRUE))
orig_auc_nxph1 <- auc(roc(model_data$Group, predict(full_model_nxph1, type="response"), quiet=TRUE))
orig_auc_combined <- auc(roc(model_data$Group, predict(full_model_combined, type="response"), quiet=TRUE))

cat(sprintf("\n%-15s %12s %12s %12s\n", "Model", "Original AUC", "LOOCV AUC", "Difference"))
cat(paste(rep("-", 55), collapse = ""), "\n")
cat(sprintf("%-15s %12.4f %12.4f %+11.4f\n", "BNC1", orig_auc_bnc1, auc(roc_bnc1_cv), auc(roc_bnc1_cv) - orig_auc_bnc1))
cat(sprintf("%-15s %12.4f %12.4f %+11.4f\n", "NXPH1", orig_auc_nxph1, auc(roc_nxph1_cv), auc(roc_nxph1_cv) - orig_auc_nxph1))
cat(sprintf("%-15s %12.4f %12.4f %+11.4f\n", "Combined", orig_auc_combined, auc(roc_combined_cv), auc(roc_combined_cv) - orig_auc_combined))

# Overfitting assessment
cat("\n=== Overfitting Assessment ===\n")
overfit_bnc1 <- (orig_auc_bnc1 - auc(roc_bnc1_cv)) / orig_auc_bnc1 * 100
overfit_nxph1 <- (orig_auc_nxph1 - auc(roc_nxph1_cv)) / orig_auc_nxph1 * 100
overfit_combined <- (orig_auc_combined - auc(roc_combined_cv)) / orig_auc_combined * 100

cat(sprintf("BNC1:     %.2f%% performance drop (", overfit_bnc1))
if(overfit_bnc1 < 5) cat("Excellent stability)\n") else if(overfit_bnc1 < 10) cat("Good stability)\n") else cat("Potential overfitting)\n")

cat(sprintf("NXPH1:    %.2f%% performance drop (", overfit_nxph1))
if(overfit_nxph1 < 5) cat("Excellent stability)\n") else if(overfit_nxph1 < 10) cat("Good stability)\n") else cat("Potential overfitting)\n")

cat(sprintf("Combined: %.2f%% performance drop (", overfit_combined))
if(overfit_combined < 5) cat("Excellent stability)\n") else if(overfit_combined < 10) cat("Good stability)\n") else cat("Potential overfitting)\n")

cat("\n=== Interpretation ===\n")
cat("- LOOCV AUC close to Original AUC = Model is stable and generalizable\n")
cat("- Large drop (>10%) = Potential overfitting, model may not generalize well\n")
cat("- Small drop (<5%) = Excellent stability, model is robust\n")

# Save results
loocv_summary <- data.frame(
  Model = c("BNC1", "NXPH1", "Combined"),
  Original_AUC = c(orig_auc_bnc1, orig_auc_nxph1, orig_auc_combined),
  LOOCV_AUC = c(auc(roc_bnc1_cv), auc(roc_nxph1_cv), auc(roc_combined_cv)),
  Accuracy = c(metrics_bnc1$Accuracy, metrics_nxph1$Accuracy, metrics_combined$Accuracy),
  Sensitivity = c(metrics_bnc1$Sensitivity, metrics_nxph1$Sensitivity, metrics_combined$Sensitivity),
  Specificity = c(metrics_bnc1$Specificity, metrics_nxph1$Specificity, metrics_combined$Specificity)
)
fwrite(loocv_summary, file.path(DMRSEQ_DIR, "loocv_results.tsv"), sep = "\t")

cat("\nResults saved to: loocv_results.tsv\n")
cat("\n=== LOOCV Analysis Complete ===\n")
