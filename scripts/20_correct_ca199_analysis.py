#!/usr/bin/env python3
"""
CA19-9 Contribution Analysis - Corrected with Control CA19-9 Data
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, roc_curve, accuracy_score, 
                             recall_score, precision_score, confusion_matrix)
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "ca199_corrected_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300

print("=" * 80)
print("CA19-9 Contribution Analysis - WITH CONTROL CA19-9 DATA")
print("=" * 80)

# 1. Load clinical data with CA19-9 for both groups
print("\n[1] Loading clinical data with CA19-9...")

patients_clinical = pd.read_excel(BASE_DIR / "data" / "환자군과 대조군260117.xlsx", sheet_name='환자군')
controls_clinical = pd.read_excel(BASE_DIR / "data" / "환자군과 대조군260117.xlsx", sheet_name='대조군')

print(f"Patients: {len(patients_clinical)} (CA19-9 available: {(~patients_clinical['CA 19-9'].isna()).sum()})")
print(f"Controls: {len(controls_clinical)} (CA19-9 available: {(~controls_clinical['CA 19-9'].isna()).sum()})")

print(f"\nCA19-9 Statistics:")
print(f"  Patients: mean={patients_clinical['CA 19-9'].mean():.1f}, median={patients_clinical['CA 19-9'].median():.1f}")
print(f"  Controls: mean={controls_clinical['CA 19-9'].mean():.1f}, median={controls_clinical['CA 19-9'].median():.1f}")
print(f"  Difference: {patients_clinical['CA 19-9'].mean() / controls_clinical['CA 19-9'].mean():.1f}x higher in patients")

# Statistical test
stat, pval = stats.mannwhitneyu(
    patients_clinical['CA 19-9'].dropna(), 
    controls_clinical['CA 19-9'].dropna(), 
    alternative='two-sided'
)
print(f"\nMann-Whitney U test: p-value = {pval:.2e}")
print(f"  -> {'Highly significant difference!' if pval < 0.001 else 'Significant difference' if pval < 0.05 else 'No significant difference'}")

# 2. Load methylation data
print("\n[2] Loading methylation data...")
integrated_data = pd.read_csv(RESULTS_DIR / "clinical_prediction" / "integrated_data.csv")

# Merge with correct CA19-9 data
patients_clinical['sample_id'] = patients_clinical['환자 번호']
controls_clinical['sample_id'] = controls_clinical['환자번호']

# Update CA19-9 in integrated data
integrated_corrected = integrated_data.copy()

# Update patient CA19-9
for idx, row in patients_clinical.iterrows():
    mask = (integrated_corrected['sample_id'] == row['sample_id']) & (integrated_corrected['label'] == 1)
    if mask.any():
        integrated_corrected.loc[mask, 'ca19_9'] = row['CA 19-9']

# Update control CA19-9
for idx, row in controls_clinical.iterrows():
    mask = (integrated_corrected['sample_id'] == row['sample_id']) & (integrated_corrected['label'] == 0)
    if mask.any():
        integrated_corrected.loc[mask, 'ca19_9'] = row['CA 19-9']

print(f"\nCorrected integrated data:")
print(f"  Total samples: {len(integrated_corrected)}")
print(f"  CA19-9 missing: {integrated_corrected['ca19_9'].isna().sum()}")
print(f"  CA19-9 available: {(~integrated_corrected['ca19_9'].isna()).sum()}")

# Save corrected data
integrated_corrected.to_csv(OUTPUT_DIR / "integrated_data_with_ca199.csv", index=False)
print(f"\nSaved corrected data: {OUTPUT_DIR / 'integrated_data_with_ca199.csv'}")

# 3. Model comparison
print("\n[3] Model Comparison with Correct CA19-9 Data...")
print("=" * 80)

dmr_cols = [col for col in integrated_corrected.columns if col.startswith('DMR_')]
clinical_cols = ['age', 'sex', 'ca19_9']

X_full = integrated_corrected[dmr_cols + clinical_cols].copy()
y = integrated_corrected['label']

# Handle missing values
for col in clinical_cols:
    if col == 'sex':
        X_full[col] = X_full[col].fillna(X_full[col].mode()[0] if len(X_full[col].mode()) > 0 else 1)
    else:
        X_full[col] = X_full[col].fillna(X_full[col].median())

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Model 1: Methylation only
print("\nModel 1: Methylation Only")
X_meth = X_full[dmr_cols].copy()
scaler_meth = StandardScaler()
X_meth_scaled = scaler_meth.fit_transform(X_meth)

lr_meth = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_meth = cross_val_predict(lr_meth, X_meth_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_meth = cross_val_predict(lr_meth, X_meth_scaled, y, cv=cv)

acc_meth = accuracy_score(y, y_pred_meth)
auc_meth = roc_auc_score(y, y_prob_meth)
sens_meth = recall_score(y, y_pred_meth)
spec_meth = recall_score(y, y_pred_meth, pos_label=0)

print(f"  Accuracy: {acc_meth:.3f}, AUC: {auc_meth:.3f}, Sensitivity: {sens_meth:.3f}, Specificity: {spec_meth:.3f}")

# Model 2: Methylation + Age + Sex
print("\nModel 2: Methylation + Age + Sex")
X_no_ca199 = X_full[dmr_cols + ['age', 'sex']].copy()
scaler_no_ca199 = StandardScaler()
X_no_ca199_scaled = scaler_no_ca199.fit_transform(X_no_ca199)

lr_no_ca199 = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_no_ca199 = cross_val_predict(lr_no_ca199, X_no_ca199_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_no_ca199 = cross_val_predict(lr_no_ca199, X_no_ca199_scaled, y, cv=cv)

acc_no_ca199 = accuracy_score(y, y_pred_no_ca199)
auc_no_ca199 = roc_auc_score(y, y_prob_no_ca199)
sens_no_ca199 = recall_score(y, y_pred_no_ca199)
spec_no_ca199 = recall_score(y, y_pred_no_ca199, pos_label=0)

print(f"  Accuracy: {acc_no_ca199:.3f}, AUC: {auc_no_ca199:.3f}, Sensitivity: {sens_no_ca199:.3f}, Specificity: {spec_no_ca199:.3f}")

# Model 3: Methylation + Age + Sex + CA19-9
print("\nModel 3: Methylation + Age + Sex + CA19-9")
X_with_ca199 = X_full[dmr_cols + clinical_cols].copy()
scaler_with_ca199 = StandardScaler()
X_with_ca199_scaled = scaler_with_ca199.fit_transform(X_with_ca199)

lr_with_ca199 = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_with_ca199 = cross_val_predict(lr_with_ca199, X_with_ca199_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_with_ca199 = cross_val_predict(lr_with_ca199, X_with_ca199_scaled, y, cv=cv)

acc_with_ca199 = accuracy_score(y, y_pred_with_ca199)
auc_with_ca199 = roc_auc_score(y, y_prob_with_ca199)
sens_with_ca199 = recall_score(y, y_pred_with_ca199)
spec_with_ca199 = recall_score(y, y_pred_with_ca199, pos_label=0)

print(f"  Accuracy: {acc_with_ca199:.3f}, AUC: {auc_with_ca199:.3f}, Sensitivity: {sens_with_ca199:.3f}, Specificity: {spec_with_ca199:.3f}")

# Model 4: CA19-9 only (for comparison)
print("\nModel 4: CA19-9 Only")
X_ca199_only = X_full[['ca19_9']].copy()
scaler_ca199_only = StandardScaler()
X_ca199_only_scaled = scaler_ca199_only.fit_transform(X_ca199_only)

lr_ca199_only = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_ca199_only = cross_val_predict(lr_ca199_only, X_ca199_only_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_ca199_only = cross_val_predict(lr_ca199_only, X_ca199_only_scaled, y, cv=cv)

acc_ca199_only = accuracy_score(y, y_pred_ca199_only)
auc_ca199_only = roc_auc_score(y, y_prob_ca199_only)
sens_ca199_only = recall_score(y, y_pred_ca199_only)
spec_ca199_only = recall_score(y, y_pred_ca199_only, pos_label=0)

print(f"  Accuracy: {acc_ca199_only:.3f}, AUC: {auc_ca199_only:.3f}, Sensitivity: {sens_ca199_only:.3f}, Specificity: {spec_ca199_only:.3f}")

# 4. Analyze CA19-9 contribution
print("\n" + "=" * 80)
print("CA19-9 CONTRIBUTION ANALYSIS")
print("=" * 80)

print(f"\nPerformance Improvements:")
print(f"  Adding CA19-9 to (Methylation + Age + Sex):")
print(f"    AUC: {auc_no_ca199:.3f} -> {auc_with_ca199:.3f} (Δ = +{auc_with_ca199-auc_no_ca199:.3f})")
print(f"    Accuracy: {acc_no_ca199*100:.1f}% -> {acc_with_ca199*100:.1f}% (Δ = +{(acc_with_ca199-acc_no_ca199)*100:.1f}%)")
print(f"    Sensitivity: {sens_no_ca199*100:.1f}% -> {sens_with_ca199*100:.1f}% (Δ = +{(sens_with_ca199-sens_no_ca199)*100:.1f}%)")
print(f"    Specificity: {spec_no_ca199*100:.1f}% -> {spec_with_ca199*100:.1f}% (Δ = +{(spec_with_ca199-spec_no_ca199)*100:.1f}%)")

# 5. Identify cases where CA19-9 helps
print("\n[5] Cases where CA19-9 changes prediction...")

results_df = pd.DataFrame({
    'sample_id': integrated_corrected['sample_id'],
    'true_label': y,
    'group': integrated_corrected['group'],
    'stage': integrated_corrected['stage'],
    'ca19_9': X_full['ca19_9'],
    'prob_no_ca199': y_prob_no_ca199,
    'prob_with_ca199': y_prob_with_ca199,
    'pred_no_ca199': y_pred_no_ca199,
    'pred_with_ca199': y_pred_with_ca199
})

results_df['prediction_changed'] = results_df['pred_no_ca199'] != results_df['pred_with_ca199']
results_df['prob_change'] = results_df['prob_with_ca199'] - results_df['prob_no_ca199']
results_df['correct_no_ca199'] = results_df['pred_no_ca199'] == results_df['true_label']
results_df['correct_with_ca199'] = results_df['pred_with_ca199'] == results_df['true_label']
results_df['ca199_helped'] = (~results_df['correct_no_ca199']) & results_df['correct_with_ca199']
results_df['ca199_hurt'] = results_df['correct_no_ca199'] & (~results_df['correct_with_ca199'])

n_changed = results_df['prediction_changed'].sum()
n_helped = results_df['ca199_helped'].sum()
n_hurt = results_df['ca199_hurt'].sum()

print(f"\nPrediction changes:")
print(f"  Total predictions changed: {n_changed}")
print(f"  CA19-9 helped (wrong -> correct): {n_helped}")
print(f"  CA19-9 hurt (correct -> wrong): {n_hurt}")
print(f"  Net benefit: {n_helped - n_hurt}")

if n_helped > 0:
    print(f"\nCases where CA19-9 HELPED:")
    helped = results_df[results_df['ca199_helped']].sort_values('prob_change', ascending=False)
    for idx, row in helped.iterrows():
        label = "Cancer" if row['true_label'] == 1 else "Control"
        stage_str = f"Stage {int(row['stage'])}" if pd.notna(row['stage']) else "N/A"
        print(f"  Sample {int(row['sample_id'])}: {label} ({stage_str})")
        print(f"    CA19-9: {row['ca19_9']:.1f}")
        print(f"    Prob: {row['prob_no_ca199']:.3f} -> {row['prob_with_ca199']:.3f} (Δ={row['prob_change']:+.3f})")
        print(f"    Prediction: {row['pred_no_ca199']} -> {row['pred_with_ca199']}")

if n_hurt > 0:
    print(f"\nCases where CA19-9 HURT:")
    hurt = results_df[results_df['ca199_hurt']].sort_values('prob_change')
    for idx, row in hurt.iterrows():
        label = "Cancer" if row['true_label'] == 1 else "Control"
        stage_str = f"Stage {int(row['stage'])}" if pd.notna(row['stage']) else "N/A"
        print(f"  Sample {int(row['sample_id'])}: {label} ({stage_str})")
        print(f"    CA19-9: {row['ca19_9']:.1f}")
        print(f"    Prob: {row['prob_no_ca199']:.3f} -> {row['prob_with_ca199']:.3f} (Δ={row['prob_change']:+.3f})")
        print(f"    Prediction: {row['pred_no_ca199']} -> {row['pred_with_ca199']}")

results_df.to_csv(OUTPUT_DIR / "ca199_contribution_details.csv", index=False)

# 6. Visualizations
print("\n[6] Generating visualizations...")

# 6-1. ROC Curves
fig, ax = plt.subplots(figsize=(10, 8))

fpr_meth, tpr_meth, _ = roc_curve(y, y_prob_meth)
fpr_ca199, tpr_ca199, _ = roc_curve(y, y_prob_ca199_only)
fpr_no_ca199, tpr_no_ca199, _ = roc_curve(y, y_prob_no_ca199)
fpr_with_ca199, tpr_with_ca199, _ = roc_curve(y, y_prob_with_ca199)

ax.plot(fpr_ca199, tpr_ca199, label=f'CA19-9 only (AUC={auc_ca199_only:.3f})', 
        linewidth=2, color='#9b59b6', linestyle='--')
ax.plot(fpr_meth, tpr_meth, label=f'Methylation only (AUC={auc_meth:.3f})', 
        linewidth=2.5, color='#3498db')
ax.plot(fpr_no_ca199, tpr_no_ca199, label=f'Methylation + Age + Sex (AUC={auc_no_ca199:.3f})', 
        linewidth=2.5, color='#e67e22')
ax.plot(fpr_with_ca199, tpr_with_ca199, label=f'Methylation + Age + Sex + CA19-9 (AUC={auc_with_ca199:.3f})', 
        linewidth=3, color='#2ecc71')
ax.plot([0, 1], [0, 1], 'k--', label='Random', linewidth=2)

ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
ax.set_title('ROC Curves: CA19-9 Contribution to Diagnosis', fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=11, loc='lower right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curves_ca199_corrected.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'roc_curves_ca199_corrected.png'}")
plt.close()

# 6-2. CA19-9 Distribution
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

cancer_ca199 = integrated_corrected[integrated_corrected['label']==1]['ca19_9'].dropna()
control_ca199 = integrated_corrected[integrated_corrected['label']==0]['ca19_9'].dropna()

# Box plot
axes[0].boxplot([control_ca199, cancer_ca199], labels=['Control', 'Cancer'])
axes[0].set_ylabel('CA19-9 Level', fontsize=12, fontweight='bold')
axes[0].set_title(f'CA19-9 Distribution\n(Control: {control_ca199.mean():.1f}, Cancer: {cancer_ca199.mean():.1f})', 
                  fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3, axis='y')

# Histogram
axes[1].hist(control_ca199, bins=20, alpha=0.6, label=f'Control (n={len(control_ca199)})', 
             color='blue', edgecolor='black')
axes[1].hist(cancer_ca199, bins=30, alpha=0.6, label=f'Cancer (n={len(cancer_ca199)})', 
             color='red', edgecolor='black')
axes[1].set_xlabel('CA19-9 Level', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Frequency', fontsize=12, fontweight='bold')
axes[1].set_title('CA19-9 Distribution Histogram', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=11)
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "ca199_distribution_corrected.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'ca199_distribution_corrected.png'}")
plt.close()

# 7. Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

summary = pd.DataFrame({
    'Model': ['CA19-9 only', 'Methylation only', 'Methylation + Age + Sex', 
              'Methylation + Age + Sex + CA19-9'],
    'Accuracy': [acc_ca199_only, acc_meth, acc_no_ca199, acc_with_ca199],
    'AUC': [auc_ca199_only, auc_meth, auc_no_ca199, auc_with_ca199],
    'Sensitivity': [sens_ca199_only, sens_meth, sens_no_ca199, sens_with_ca199],
    'Specificity': [spec_ca199_only, spec_meth, spec_no_ca199, spec_with_ca199]
})

summary.to_csv(OUTPUT_DIR / "ca199_contribution_summary.csv", index=False)
print(f"\nSummary saved: {OUTPUT_DIR / 'ca199_contribution_summary.csv'}")
print("\n" + summary.to_string(index=False))

print(f"\nKEY FINDINGS:")
print(f"1. CA19-9 alone: AUC = {auc_ca199_only:.3f}")
print(f"2. Methylation alone: AUC = {auc_meth:.3f}")
print(f"3. Adding CA19-9 to full model improved AUC by {auc_with_ca199-auc_no_ca199:.3f}")
print(f"4. CA19-9 helped {n_helped} cases, hurt {n_hurt} cases (net: {n_helped-n_hurt})")
print(f"5. CA19-9 is {cancer_ca199.mean()/control_ca199.mean():.1f}x higher in cancer vs control (p={pval:.2e})")

print(f"\nResults saved to: {OUTPUT_DIR}")
print("=" * 80)
