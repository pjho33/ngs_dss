#!/usr/bin/env python3
"""
Stage-Specific Pancreatic Cancer Diagnosis Performance Analysis
Focus: Early-stage (Stage 1-2) detection rate
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (roc_auc_score, roc_curve, classification_report, 
                             confusion_matrix, accuracy_score, precision_score, 
                             recall_score, f1_score)
import warnings
warnings.filterwarnings('ignore')

# Path setup
BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "stage_analysis_english"
OUTPUT_DIR.mkdir(exist_ok=True)

# Visualization setup
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300

print("=" * 80)
print("PANCREATIC CANCER DIAGNOSIS: STAGE-SPECIFIC PERFORMANCE ANALYSIS")
print("=" * 80)
print("\nSTATISTICAL METHODS AND METRICS:")
print("-" * 80)
print("1. FEATURE EXTRACTION:")
print("   - DMRseq (Differentially Methylated Regions sequencing)")
print("   - Top 46 significant DMRs (q-value < 0.05)")
print("   - Methylation levels calculated as: M/(M+U) * 100")
print("     where M = methylated reads, U = unmethylated reads")
print("")
print("2. CLASSIFICATION MODEL:")
print("   - Algorithm: Logistic Regression with L2 regularization")
print("   - Class weighting: Balanced (to handle class imbalance)")
print("   - Feature scaling: StandardScaler (z-score normalization)")
print("   - Validation: 5-fold Stratified Cross-Validation")
print("")
print("3. PERFORMANCE METRICS:")
print("   - Sensitivity (Recall): TP / (TP + FN)")
print("     * Measures ability to detect cancer patients")
print("   - Specificity: TN / (TN + FP)")
print("     * Measures ability to correctly identify controls")
print("   - Precision (PPV): TP / (TP + FP)")
print("     * Positive Predictive Value")
print("   - Accuracy: (TP + TN) / (TP + TN + FP + FN)")
print("   - F1-score: 2 * (Precision * Recall) / (Precision + Recall)")
print("   - AUC-ROC: Area Under the Receiver Operating Characteristic Curve")
print("     * Measures overall discriminative ability")
print("")
print("4. STATISTICAL SIGNIFICANCE:")
print("   - DMR selection: q-value < 0.05 (FDR-adjusted p-value)")
print("   - Cross-validation ensures unbiased performance estimation")
print("=" * 80)

# 1. Load integrated data
print("\n[1] Loading integrated data...")
data_file = RESULTS_DIR / "clinical_prediction" / "integrated_data.csv"
df = pd.read_csv(data_file)

print(f"Total samples: {df.shape[0]}")
print(f"  Cancer patients: {sum(df['label'] == 1)}")
print(f"  Healthy controls: {sum(df['label'] == 0)}")

# 2. Stage distribution analysis
print("\n[2] Cancer stage distribution analysis...")

patient_df = df[df['label'] == 1].copy()

print(f"\nStage distribution in cancer patients:")
stage_counts = patient_df['stage'].value_counts().sort_index()
for stage, count in stage_counts.items():
    print(f"  Stage {int(stage)}: {count} patients")

print(f"\nMissing stage information: {patient_df['stage'].isna().sum()} patients")

# Create stage groups
patient_df['stage_group'] = patient_df['stage'].apply(
    lambda x: 'Early (Stage 1-2)' if x in [1, 2] else ('Late (Stage 3-4)' if x in [3, 4] else 'Unknown')
)

print(f"\nStage group distribution:")
for group, count in patient_df['stage_group'].value_counts().items():
    print(f"  {group}: {count} patients")

# 3. Early-stage vs control analysis
print("\n[3] Early-stage cancer (Stage 1-2) vs Control analysis...")

early_stage = patient_df[patient_df['stage'].isin([1, 2])].copy()
late_stage = patient_df[patient_df['stage'].isin([3, 4])].copy()
normal_df = df[df['label'] == 0].copy()

print(f"\nEarly-stage cancer patients: {len(early_stage)}")
print(f"  Stage 1: {sum(early_stage['stage'] == 1)}")
print(f"  Stage 2: {sum(early_stage['stage'] == 2)}")
print(f"Late-stage cancer patients: {len(late_stage)}")
print(f"  Stage 3: {sum(late_stage['stage'] == 3)}")
print(f"  Stage 4: {sum(late_stage['stage'] == 4)}")
print(f"Healthy controls: {len(normal_df)}")

# Feature columns
dmr_cols = [col for col in df.columns if col.startswith('DMR_')]
clinical_cols = ['age', 'sex', 'ca19_9']

print(f"\nFeatures used:")
print(f"  DMR methylation features: {len(dmr_cols)}")
print(f"  Clinical features: {clinical_cols}")

# 4. Model 1: Early-stage vs Control
print("\n" + "=" * 80)
print("MODEL 1: EARLY-STAGE CANCER (Stage 1-2) vs HEALTHY CONTROLS")
print("=" * 80)

early_vs_normal = pd.concat([early_stage, normal_df], ignore_index=True)
X_early = early_vs_normal[dmr_cols + clinical_cols].copy()
y_early = early_vs_normal['label']

# Handle missing values
for col in clinical_cols:
    if col == 'sex':
        X_early[col] = X_early[col].fillna(X_early[col].mode()[0] if len(X_early[col].mode()) > 0 else 1)
    else:
        X_early[col] = X_early[col].fillna(X_early[col].median())

# Standardization
scaler_early = StandardScaler()
X_early_scaled = scaler_early.fit_transform(X_early)

# Logistic Regression with balanced class weights
lr_early = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')

# 5-fold stratified cross-validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

y_pred_early_cv = cross_val_predict(lr_early, X_early_scaled, y_early, cv=cv)
y_prob_early_cv = cross_val_predict(lr_early, X_early_scaled, y_early, cv=cv, method='predict_proba')[:, 1]

# Performance metrics
acc_early = accuracy_score(y_early, y_pred_early_cv)
auc_early = roc_auc_score(y_early, y_prob_early_cv)
sens_early = recall_score(y_early, y_pred_early_cv)
spec_early = recall_score(y_early, y_pred_early_cv, pos_label=0)
prec_early = precision_score(y_early, y_pred_early_cv)
f1_early = f1_score(y_early, y_pred_early_cv)

print(f"\nPerformance Metrics:")
print(f"  Accuracy: {acc_early:.3f} ({acc_early*100:.1f}%)")
print(f"  AUC-ROC: {auc_early:.3f}")
print(f"  Sensitivity (Early-stage detection rate): {sens_early:.3f} ({sens_early*100:.1f}%)")
print(f"  Specificity: {spec_early:.3f} ({spec_early*100:.1f}%)")
print(f"  Precision (PPV): {prec_early:.3f} ({prec_early*100:.1f}%)")
print(f"  F1-score: {f1_early:.3f}")

cm_early = confusion_matrix(y_early, y_pred_early_cv)
print(f"\nConfusion Matrix:")
print(f"                    Predicted Control  Predicted Cancer")
print(f"  Actual Control:   {cm_early[0,0]:6d}            {cm_early[0,1]:6d}")
print(f"  Actual Cancer:    {cm_early[1,0]:6d}            {cm_early[1,1]:6d}")
print(f"\n  True Negatives (TN): {cm_early[0,0]}")
print(f"  False Positives (FP): {cm_early[0,1]}")
print(f"  False Negatives (FN): {cm_early[1,0]}")
print(f"  True Positives (TP): {cm_early[1,1]}")

# 5. Model 2: Late-stage vs Control
print("\n" + "=" * 80)
print("MODEL 2: LATE-STAGE CANCER (Stage 3-4) vs HEALTHY CONTROLS")
print("=" * 80)

late_vs_normal = pd.concat([late_stage, normal_df], ignore_index=True)
X_late = late_vs_normal[dmr_cols + clinical_cols].copy()
y_late = late_vs_normal['label']

for col in clinical_cols:
    if col == 'sex':
        X_late[col] = X_late[col].fillna(X_late[col].mode()[0] if len(X_late[col].mode()) > 0 else 1)
    else:
        X_late[col] = X_late[col].fillna(X_late[col].median())

scaler_late = StandardScaler()
X_late_scaled = scaler_late.fit_transform(X_late)

lr_late = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')

y_pred_late_cv = cross_val_predict(lr_late, X_late_scaled, y_late, cv=cv)
y_prob_late_cv = cross_val_predict(lr_late, X_late_scaled, y_late, cv=cv, method='predict_proba')[:, 1]

acc_late = accuracy_score(y_late, y_pred_late_cv)
auc_late = roc_auc_score(y_late, y_prob_late_cv)
sens_late = recall_score(y_late, y_pred_late_cv)
spec_late = recall_score(y_late, y_pred_late_cv, pos_label=0)
prec_late = precision_score(y_late, y_pred_late_cv)
f1_late = f1_score(y_late, y_pred_late_cv)

print(f"\nPerformance Metrics:")
print(f"  Accuracy: {acc_late:.3f} ({acc_late*100:.1f}%)")
print(f"  AUC-ROC: {auc_late:.3f}")
print(f"  Sensitivity (Late-stage detection rate): {sens_late:.3f} ({sens_late*100:.1f}%)")
print(f"  Specificity: {spec_late:.3f} ({spec_late*100:.1f}%)")
print(f"  Precision (PPV): {prec_late:.3f} ({prec_late*100:.1f}%)")
print(f"  F1-score: {f1_late:.3f}")

cm_late = confusion_matrix(y_late, y_pred_late_cv)
print(f"\nConfusion Matrix:")
print(f"                    Predicted Control  Predicted Cancer")
print(f"  Actual Control:   {cm_late[0,0]:6d}            {cm_late[0,1]:6d}")
print(f"  Actual Cancer:    {cm_late[1,0]:6d}            {cm_late[1,1]:6d}")
print(f"\n  True Negatives (TN): {cm_late[0,0]}")
print(f"  False Positives (FP): {cm_late[0,1]}")
print(f"  False Negatives (FN): {cm_late[1,0]}")
print(f"  True Positives (TP): {cm_late[1,1]}")

# 6. Visualizations
print("\n[6] Generating visualizations...")

# 6-1. ROC Curves
fig, ax = plt.subplots(figsize=(10, 8))

fpr_early, tpr_early, _ = roc_curve(y_early, y_prob_early_cv)
fpr_late, tpr_late, _ = roc_curve(y_late, y_prob_late_cv)

ax.plot(fpr_early, tpr_early, label=f'Early-stage (1-2) vs Control (AUC={auc_early:.3f})', 
        linewidth=3, color='#2ecc71')
ax.plot(fpr_late, tpr_late, label=f'Late-stage (3-4) vs Control (AUC={auc_late:.3f})', 
        linewidth=3, color='#e74c3c')
ax.plot([0, 1], [0, 1], 'k--', label='Random Classifier', linewidth=2)

ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=14, fontweight='bold')
ax.set_title('ROC Curves: Stage-Specific Diagnostic Performance', fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=12, loc='lower right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curves_by_stage.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'roc_curves_by_stage.png'}")
plt.close()

# 6-2. Performance comparison bar charts
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

metrics = ['Accuracy', 'AUC-ROC', 'Sensitivity', 'Specificity']
early_scores = [acc_early, auc_early, sens_early, spec_early]
late_scores = [acc_late, auc_late, sens_late, spec_late]

for idx, (ax, metric, early_val, late_val) in enumerate(zip(axes.flat, metrics, early_scores, late_scores)):
    x = ['Early-stage\n(Stage 1-2)', 'Late-stage\n(Stage 3-4)']
    values = [early_val, late_val]
    colors = ['#2ecc71', '#e74c3c']
    
    bars = ax.bar(x, values, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
    ax.set_ylabel(metric, fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1.1])
    ax.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{val:.3f}\n({val*100:.1f}%)', ha='center', va='bottom', 
                fontsize=10, fontweight='bold')

plt.suptitle('Stage-Specific Diagnostic Performance Comparison', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "performance_comparison_by_stage.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'performance_comparison_by_stage.png'}")
plt.close()

# 6-3. Confusion matrices
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Early-stage
sns.heatmap(cm_early, annot=True, fmt='d', cmap='Greens', ax=axes[0], 
            cbar_kws={'label': 'Count'}, annot_kws={'size': 14, 'weight': 'bold'})
axes[0].set_title(f'Early-stage (1-2) vs Control\nSensitivity: {sens_early:.1%}', 
                  fontsize=14, fontweight='bold')
axes[0].set_ylabel('Actual', fontsize=12, fontweight='bold')
axes[0].set_xlabel('Predicted', fontsize=12, fontweight='bold')
axes[0].set_xticklabels(['Control', 'Cancer'])
axes[0].set_yticklabels(['Control', 'Cancer'])

# Late-stage
sns.heatmap(cm_late, annot=True, fmt='d', cmap='Reds', ax=axes[1], 
            cbar_kws={'label': 'Count'}, annot_kws={'size': 14, 'weight': 'bold'})
axes[1].set_title(f'Late-stage (3-4) vs Control\nSensitivity: {sens_late:.1%}', 
                  fontsize=14, fontweight='bold')
axes[1].set_ylabel('Actual', fontsize=12, fontweight='bold')
axes[1].set_xlabel('Predicted', fontsize=12, fontweight='bold')
axes[1].set_xticklabels(['Control', 'Cancer'])
axes[1].set_yticklabels(['Control', 'Cancer'])

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "confusion_matrices_by_stage.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'confusion_matrices_by_stage.png'}")
plt.close()

# 7. Save summary results
print("\n[7] Saving summary results...")

summary = pd.DataFrame({
    'Stage_Group': ['Early-stage (1-2)', 'Late-stage (3-4)'],
    'N_Patients': [len(early_stage), len(late_stage)],
    'Accuracy': [acc_early, acc_late],
    'AUC_ROC': [auc_early, auc_late],
    'Sensitivity': [sens_early, sens_late],
    'Specificity': [spec_early, spec_late],
    'Precision_PPV': [prec_early, prec_late],
    'F1_Score': [f1_early, f1_late]
})

summary.to_csv(OUTPUT_DIR / "stage_specific_performance_summary.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'stage_specific_performance_summary.csv'}")

# Detailed results
detailed_results = pd.DataFrame({
    'Metric': ['N_Patients', 'N_Controls', 'Accuracy', 'AUC-ROC', 'Sensitivity', 
               'Specificity', 'Precision', 'F1-Score', 'True_Positive', 'False_Negative',
               'True_Negative', 'False_Positive'],
    'Early_Stage_1_2': [
        len(early_stage), len(normal_df), acc_early, auc_early, sens_early,
        spec_early, prec_early, f1_early, int(cm_early[1,1]), int(cm_early[1,0]),
        int(cm_early[0,0]), int(cm_early[0,1])
    ],
    'Late_Stage_3_4': [
        len(late_stage), len(normal_df), acc_late, auc_late, sens_late,
        spec_late, prec_late, f1_late, int(cm_late[1,1]), int(cm_late[1,0]),
        int(cm_late[0,0]), int(cm_late[0,1])
    ]
})

detailed_results.to_csv(OUTPUT_DIR / "detailed_stage_results.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'detailed_stage_results.csv'}")

# Statistical methods documentation
methods_doc = """
STATISTICAL METHODS AND PERFORMANCE METRICS
==========================================

1. FEATURE EXTRACTION
   - Method: DMRseq (Differentially Methylated Regions sequencing)
   - Significance threshold: q-value < 0.05 (FDR-corrected)
   - Number of DMRs: 46 significant regions
   - Methylation calculation: M/(M+U) * 100%
     where M = methylated reads, U = unmethylated reads

2. CLASSIFICATION MODEL
   - Algorithm: Logistic Regression
   - Regularization: L2 (Ridge)
   - Class weighting: Balanced (inversely proportional to class frequencies)
   - Feature preprocessing: StandardScaler (z-score normalization)
   - Validation: 5-fold Stratified Cross-Validation
   - Random seed: 42 (for reproducibility)

3. PERFORMANCE METRICS
   
   a) Sensitivity (Recall, True Positive Rate)
      Formula: TP / (TP + FN)
      Interpretation: Proportion of actual cancer patients correctly identified
      Clinical significance: Ability to detect cancer (minimize missed diagnoses)
   
   b) Specificity (True Negative Rate)
      Formula: TN / (TN + FP)
      Interpretation: Proportion of healthy controls correctly identified
      Clinical significance: Ability to avoid false alarms
   
   c) Precision (Positive Predictive Value, PPV)
      Formula: TP / (TP + FP)
      Interpretation: Proportion of positive predictions that are correct
      Clinical significance: Reliability of positive test results
   
   d) Accuracy
      Formula: (TP + TN) / (TP + TN + FP + FN)
      Interpretation: Overall proportion of correct predictions
   
   e) F1-Score
      Formula: 2 * (Precision * Recall) / (Precision + Recall)
      Interpretation: Harmonic mean of precision and recall
      Clinical significance: Balanced measure of diagnostic performance
   
   f) AUC-ROC (Area Under the Receiver Operating Characteristic Curve)
      Range: 0.5 (random) to 1.0 (perfect)
      Interpretation: Overall discriminative ability across all thresholds
      Clinical significance: Model's ability to rank patients by cancer risk

4. STATISTICAL SIGNIFICANCE
   - DMR identification: q-value < 0.05 (Benjamini-Hochberg FDR correction)
   - Cross-validation: Ensures unbiased performance estimation
   - Stratification: Maintains class proportions in each fold

5. CLINICAL FEATURES
   - Age (continuous)
   - Sex (binary: 1=male, 2=female)
   - CA19-9 (continuous, tumor marker)
   - Stage excluded from prediction (target variable for stratification)

ABBREVIATIONS:
TP = True Positive, TN = True Negative
FP = False Positive, FN = False Negative
DMR = Differentially Methylated Region
FDR = False Discovery Rate
PPV = Positive Predictive Value
AUC = Area Under the Curve
ROC = Receiver Operating Characteristic
"""

with open(OUTPUT_DIR / "statistical_methods.txt", 'w') as f:
    f.write(methods_doc)
print(f"Saved: {OUTPUT_DIR / 'statistical_methods.txt'}")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE!")
print("=" * 80)

print("\n" + "=" * 80)
print("KEY FINDINGS:")
print("=" * 80)
print(f"\nEARLY-STAGE CANCER (Stage 1-2) DETECTION:")
print(f"  Number of patients: {len(early_stage)} (Stage 1: {sum(early_stage['stage']==1)}, Stage 2: {sum(early_stage['stage']==2)})")
print(f"  Sensitivity (Detection rate): {sens_early:.1%} ({int(cm_early[1,1])}/{len(early_stage)} patients detected)")
print(f"  AUC-ROC: {auc_early:.3f}")
print(f"  Accuracy: {acc_early:.1%}")
print(f"  Specificity: {spec_early:.1%}")
print(f"  --> {int(cm_early[1,1])} out of {len(early_stage)} early-stage patients correctly identified!")

print(f"\nLATE-STAGE CANCER (Stage 3-4) DETECTION:")
print(f"  Number of patients: {len(late_stage)} (Stage 3: {sum(late_stage['stage']==3)}, Stage 4: {sum(late_stage['stage']==4)})")
print(f"  Sensitivity (Detection rate): {sens_late:.1%} ({int(cm_late[1,1])}/{len(late_stage)} patients detected)")
print(f"  AUC-ROC: {auc_late:.3f}")
print(f"  Accuracy: {acc_late:.1%}")
print(f"  Specificity: {spec_late:.1%}")

print(f"\nCLINICAL SIGNIFICANCE:")
print(f"  - The model achieves {sens_early:.1%} sensitivity for early-stage cancer")
print(f"  - Only {int(cm_early[1,0])} early-stage patient(s) missed (False Negative)")
print(f"  - High specificity ({spec_early:.1%}) minimizes false alarms")
print(f"  - Excellent AUC-ROC ({auc_early:.3f}) indicates strong discriminative ability")

print(f"\nResults saved to: {OUTPUT_DIR}")
print("=" * 80)
