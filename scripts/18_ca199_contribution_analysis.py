#!/usr/bin/env python3
"""
CA19-9의 진단 기여도 분석
메틸화 DNA만으로 진단이 어려운 경우 CA19-9가 도움이 되는지 평가
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
import warnings
warnings.filterwarnings('ignore')

# Path setup
BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "ca199_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300

print("=" * 80)
print("CA19-9 Contribution Analysis for Pancreatic Cancer Diagnosis")
print("=" * 80)

# 1. Load integrated data
print("\n[1] Loading integrated data...")
data_file = RESULTS_DIR / "clinical_prediction" / "integrated_data.csv"
df = pd.read_csv(data_file)

print(f"Total samples: {df.shape[0]}")
print(f"  Cancer patients: {sum(df['label'] == 1)}")
print(f"  Healthy controls: {sum(df['label'] == 0)}")

# Feature columns
dmr_cols = [col for col in df.columns if col.startswith('DMR_')]
clinical_cols = ['age', 'sex', 'ca19_9']

print(f"\nFeatures:")
print(f"  DMR methylation features: {len(dmr_cols)}")
print(f"  Clinical features: {clinical_cols}")

# 2. Prepare data
print("\n[2] Preparing data...")

# Handle missing values
X_full = df[dmr_cols + clinical_cols].copy()
y = df['label']

for col in clinical_cols:
    if col == 'sex':
        X_full[col] = X_full[col].fillna(X_full[col].mode()[0] if len(X_full[col].mode()) > 0 else 1)
    else:
        X_full[col] = X_full[col].fillna(X_full[col].median())

# CA19-9 statistics
print(f"\nCA19-9 Statistics:")
print(f"  Cancer patients (mean ± std): {X_full[y==1]['ca19_9'].mean():.1f} ± {X_full[y==1]['ca19_9'].std():.1f}")
print(f"  Controls (mean ± std): {X_full[y==0]['ca19_9'].mean():.1f} ± {X_full[y==0]['ca19_9'].std():.1f}")
print(f"  Missing values: {df['ca19_9'].isna().sum()}")

# 3. Model Comparison
print("\n[3] Model Comparison...")
print("=" * 80)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Model 1: Methylation only
print("\nModel 1: Methylation Only (46 DMRs)")
print("-" * 80)

X_meth = X_full[dmr_cols].copy()
scaler_meth = StandardScaler()
X_meth_scaled = scaler_meth.fit_transform(X_meth)

lr_meth = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_pred_meth = cross_val_predict(lr_meth, X_meth_scaled, y, cv=cv)
y_prob_meth = cross_val_predict(lr_meth, X_meth_scaled, y, cv=cv, method='predict_proba')[:, 1]

acc_meth = accuracy_score(y, y_pred_meth)
auc_meth = roc_auc_score(y, y_prob_meth)
sens_meth = recall_score(y, y_pred_meth)
spec_meth = recall_score(y, y_pred_meth, pos_label=0)
prec_meth = precision_score(y, y_pred_meth)

print(f"Accuracy: {acc_meth:.3f} ({acc_meth*100:.1f}%)")
print(f"AUC-ROC: {auc_meth:.3f}")
print(f"Sensitivity: {sens_meth:.3f} ({sens_meth*100:.1f}%)")
print(f"Specificity: {spec_meth:.3f} ({spec_meth*100:.1f}%)")
print(f"Precision: {prec_meth:.3f} ({prec_meth*100:.1f}%)")

# Model 2: Methylation + Age + Sex (without CA19-9)
print("\nModel 2: Methylation + Age + Sex (without CA19-9)")
print("-" * 80)

X_no_ca199 = X_full[dmr_cols + ['age', 'sex']].copy()
scaler_no_ca199 = StandardScaler()
X_no_ca199_scaled = scaler_no_ca199.fit_transform(X_no_ca199)

lr_no_ca199 = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_pred_no_ca199 = cross_val_predict(lr_no_ca199, X_no_ca199_scaled, y, cv=cv)
y_prob_no_ca199 = cross_val_predict(lr_no_ca199, X_no_ca199_scaled, y, cv=cv, method='predict_proba')[:, 1]

acc_no_ca199 = accuracy_score(y, y_pred_no_ca199)
auc_no_ca199 = roc_auc_score(y, y_prob_no_ca199)
sens_no_ca199 = recall_score(y, y_pred_no_ca199)
spec_no_ca199 = recall_score(y, y_pred_no_ca199, pos_label=0)
prec_no_ca199 = precision_score(y, y_pred_no_ca199)

print(f"Accuracy: {acc_no_ca199:.3f} ({acc_no_ca199*100:.1f}%)")
print(f"AUC-ROC: {auc_no_ca199:.3f}")
print(f"Sensitivity: {sens_no_ca199:.3f} ({sens_no_ca199*100:.1f}%)")
print(f"Specificity: {spec_no_ca199:.3f} ({spec_no_ca199*100:.1f}%)")
print(f"Precision: {prec_no_ca199:.3f} ({prec_no_ca199*100:.1f}%)")

# Model 3: Methylation + Age + Sex + CA19-9
print("\nModel 3: Methylation + Age + Sex + CA19-9")
print("-" * 80)

X_with_ca199 = X_full[dmr_cols + clinical_cols].copy()
scaler_with_ca199 = StandardScaler()
X_with_ca199_scaled = scaler_with_ca199.fit_transform(X_with_ca199)

lr_with_ca199 = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_pred_with_ca199 = cross_val_predict(lr_with_ca199, X_with_ca199_scaled, y, cv=cv)
y_prob_with_ca199 = cross_val_predict(lr_with_ca199, X_with_ca199_scaled, y, cv=cv, method='predict_proba')[:, 1]

acc_with_ca199 = accuracy_score(y, y_pred_with_ca199)
auc_with_ca199 = roc_auc_score(y, y_prob_with_ca199)
sens_with_ca199 = recall_score(y, y_pred_with_ca199)
spec_with_ca199 = recall_score(y, y_pred_with_ca199, pos_label=0)
prec_with_ca199 = precision_score(y, y_pred_with_ca199)

print(f"Accuracy: {acc_with_ca199:.3f} ({acc_with_ca199*100:.1f}%)")
print(f"AUC-ROC: {auc_with_ca199:.3f}")
print(f"Sensitivity: {sens_with_ca199:.3f} ({sens_with_ca199*100:.1f}%)")
print(f"Specificity: {spec_with_ca199:.3f} ({spec_with_ca199*100:.1f}%)")
print(f"Precision: {prec_with_ca199:.3f} ({prec_with_ca199*100:.1f}%)")

# Calculate improvement
print("\n" + "=" * 80)
print("CA19-9 CONTRIBUTION ANALYSIS")
print("=" * 80)

auc_improvement = auc_with_ca199 - auc_no_ca199
acc_improvement = acc_with_ca199 - acc_no_ca199
sens_improvement = sens_with_ca199 - sens_no_ca199

print(f"\nPerformance Improvement by Adding CA19-9:")
print(f"  AUC-ROC: {auc_no_ca199:.3f} -> {auc_with_ca199:.3f} (Δ = +{auc_improvement:.3f})")
print(f"  Accuracy: {acc_no_ca199*100:.1f}% -> {acc_with_ca199*100:.1f}% (Δ = +{acc_improvement*100:.1f}%)")
print(f"  Sensitivity: {sens_no_ca199*100:.1f}% -> {sens_with_ca199*100:.1f}% (Δ = +{sens_improvement*100:.1f}%)")

# 4. Identify cases where CA19-9 helps
print("\n[4] Identifying cases where CA19-9 helps diagnosis...")
print("=" * 80)

# Get prediction probabilities for both models
results_df = pd.DataFrame({
    'sample_id': df['sample_id'],
    'true_label': y,
    'group': df['group'],
    'ca19_9': X_full['ca19_9'],
    'prob_without_ca199': y_prob_no_ca199,
    'prob_with_ca199': y_prob_with_ca199,
    'pred_without_ca199': y_pred_no_ca199,
    'pred_with_ca199': y_pred_with_ca199
})

# Cases where predictions differ
results_df['prediction_changed'] = results_df['pred_without_ca199'] != results_df['pred_with_ca199']
results_df['prob_change'] = results_df['prob_with_ca199'] - results_df['prob_without_ca199']

# Correctly diagnosed only with CA19-9
results_df['correct_without'] = results_df['pred_without_ca199'] == results_df['true_label']
results_df['correct_with'] = results_df['pred_with_ca199'] == results_df['true_label']
results_df['ca199_helped'] = (~results_df['correct_without']) & results_df['correct_with']
results_df['ca199_hurt'] = results_df['correct_without'] & (~results_df['correct_with'])

n_helped = results_df['ca199_helped'].sum()
n_hurt = results_df['ca199_hurt'].sum()
n_changed = results_df['prediction_changed'].sum()

print(f"\nPrediction Changes:")
print(f"  Total predictions changed: {n_changed}")
print(f"  CA19-9 helped (wrong -> correct): {n_helped}")
print(f"  CA19-9 hurt (correct -> wrong): {n_hurt}")
print(f"  Net benefit: {n_helped - n_hurt}")

# Cases where CA19-9 helped
if n_helped > 0:
    print(f"\n{n_helped} cases where CA19-9 helped diagnosis:")
    helped_cases = results_df[results_df['ca199_helped']].copy()
    helped_cases = helped_cases.sort_values('prob_change', ascending=False)
    
    for idx, row in helped_cases.iterrows():
        label_str = "Cancer" if row['true_label'] == 1 else "Control"
        print(f"  Sample {int(row['sample_id'])}: {label_str}")
        print(f"    CA19-9: {row['ca19_9']:.1f}")
        print(f"    Prob without CA19-9: {row['prob_without_ca199']:.3f} -> Pred: {row['pred_without_ca199']}")
        print(f"    Prob with CA19-9: {row['prob_with_ca199']:.3f} -> Pred: {row['pred_with_ca199']}")
        print(f"    Probability change: +{row['prob_change']:.3f}")
        print()

# Difficult cases (low confidence without CA19-9)
results_df['confidence_without'] = np.abs(results_df['prob_without_ca199'] - 0.5)
results_df['confidence_with'] = np.abs(results_df['prob_with_ca199'] - 0.5)
results_df['confidence_improved'] = results_df['confidence_with'] > results_df['confidence_without']

difficult_cases = results_df[results_df['confidence_without'] < 0.2].copy()
print(f"\nDifficult cases (low confidence without CA19-9, prob ~ 0.5):")
print(f"  Total difficult cases: {len(difficult_cases)}")
print(f"  Confidence improved with CA19-9: {difficult_cases['confidence_improved'].sum()}")

if len(difficult_cases) > 0:
    print(f"\n  Top 10 difficult cases:")
    difficult_cases = difficult_cases.sort_values('confidence_without')
    for idx, row in difficult_cases.head(10).iterrows():
        label_str = "Cancer" if row['true_label'] == 1 else "Control"
        correct_str = "✓" if row['correct_with'] else "✗"
        print(f"    Sample {int(row['sample_id'])} ({label_str}): {correct_str}")
        print(f"      CA19-9: {row['ca19_9']:.1f}")
        print(f"      Confidence without CA19-9: {row['confidence_without']:.3f}")
        print(f"      Confidence with CA19-9: {row['confidence_with']:.3f}")

# Save results
results_df.to_csv(OUTPUT_DIR / "ca199_contribution_details.csv", index=False)
print(f"\nDetailed results saved: {OUTPUT_DIR / 'ca199_contribution_details.csv'}")

# 5. Visualizations
print("\n[5] Generating visualizations...")

# 5-1. ROC Curves Comparison
fig, ax = plt.subplots(figsize=(10, 8))

fpr_meth, tpr_meth, _ = roc_curve(y, y_prob_meth)
fpr_no_ca199, tpr_no_ca199, _ = roc_curve(y, y_prob_no_ca199)
fpr_with_ca199, tpr_with_ca199, _ = roc_curve(y, y_prob_with_ca199)

ax.plot(fpr_meth, tpr_meth, label=f'Methylation only (AUC={auc_meth:.3f})', 
        linewidth=2.5, color='#3498db')
ax.plot(fpr_no_ca199, tpr_no_ca199, label=f'Methylation + Age + Sex (AUC={auc_no_ca199:.3f})', 
        linewidth=2.5, color='#e67e22')
ax.plot(fpr_with_ca199, tpr_with_ca199, label=f'Methylation + Age + Sex + CA19-9 (AUC={auc_with_ca199:.3f})', 
        linewidth=3, color='#2ecc71')
ax.plot([0, 1], [0, 1], 'k--', label='Random', linewidth=2)

ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
ax.set_title('ROC Curves: Impact of CA19-9 on Diagnostic Performance', 
             fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=11, loc='lower right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curves_ca199_comparison.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'roc_curves_ca199_comparison.png'}")
plt.close()

# 5-2. Performance Metrics Comparison
fig, ax = plt.subplots(figsize=(12, 6))

metrics = ['Accuracy', 'AUC-ROC', 'Sensitivity', 'Specificity']
meth_scores = [acc_meth, auc_meth, sens_meth, spec_meth]
no_ca199_scores = [acc_no_ca199, auc_no_ca199, sens_no_ca199, spec_no_ca199]
with_ca199_scores = [acc_with_ca199, auc_with_ca199, sens_with_ca199, spec_with_ca199]

x = np.arange(len(metrics))
width = 0.25

bars1 = ax.bar(x - width, meth_scores, width, label='Methylation only', 
               color='#3498db', alpha=0.8, edgecolor='black', linewidth=1.5)
bars2 = ax.bar(x, no_ca199_scores, width, label='+ Age + Sex', 
               color='#e67e22', alpha=0.8, edgecolor='black', linewidth=1.5)
bars3 = ax.bar(x + width, with_ca199_scores, width, label='+ Age + Sex + CA19-9', 
               color='#2ecc71', alpha=0.8, edgecolor='black', linewidth=1.5)

ax.set_ylabel('Score', fontsize=13, fontweight='bold')
ax.set_title('Performance Comparison: CA19-9 Contribution', fontsize=15, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=12)
ax.legend(fontsize=11)
ax.set_ylim([0, 1.1])
ax.grid(True, alpha=0.3, axis='y')

# Add value labels
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "performance_comparison_ca199.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'performance_comparison_ca199.png'}")
plt.close()

# 5-3. CA19-9 Distribution by Group
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Box plot
cancer_ca199 = X_full[y==1]['ca19_9']
control_ca199 = X_full[y==0]['ca19_9']

axes[0].boxplot([control_ca199, cancer_ca199], labels=['Control', 'Cancer'])
axes[0].set_ylabel('CA19-9 Level', fontsize=12, fontweight='bold')
axes[0].set_title('CA19-9 Distribution by Group', fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3, axis='y')

# Histogram
axes[1].hist(control_ca199, bins=20, alpha=0.6, label='Control', color='blue', edgecolor='black')
axes[1].hist(cancer_ca199, bins=20, alpha=0.6, label='Cancer', color='red', edgecolor='black')
axes[1].set_xlabel('CA19-9 Level', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Frequency', fontsize=12, fontweight='bold')
axes[1].set_title('CA19-9 Distribution Histogram', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=11)
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "ca199_distribution.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'ca199_distribution.png'}")
plt.close()

# 6. Summary Report
print("\n" + "=" * 80)
print("SUMMARY: CA19-9 CONTRIBUTION TO PANCREATIC CANCER DIAGNOSIS")
print("=" * 80)

summary = {
    'Model': ['Methylation only', 'Methylation + Age + Sex', 'Methylation + Age + Sex + CA19-9'],
    'Accuracy': [acc_meth, acc_no_ca199, acc_with_ca199],
    'AUC_ROC': [auc_meth, auc_no_ca199, auc_with_ca199],
    'Sensitivity': [sens_meth, sens_no_ca199, sens_with_ca199],
    'Specificity': [spec_meth, spec_no_ca199, spec_with_ca199],
    'Precision': [prec_meth, prec_no_ca199, prec_with_ca199]
}

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUTPUT_DIR / "ca199_contribution_summary.csv", index=False)
print(f"\nSummary saved: {OUTPUT_DIR / 'ca199_contribution_summary.csv'}")

print("\n" + summary_df.to_string(index=False))

print(f"\nKEY FINDINGS:")
print(f"1. CA19-9 improved AUC-ROC by {auc_improvement:.3f} ({auc_improvement/auc_no_ca199*100:.1f}%)")
print(f"2. CA19-9 helped correctly diagnose {n_helped} additional cases")
print(f"3. CA19-9 improved confidence in {difficult_cases['confidence_improved'].sum()} difficult cases")
print(f"4. Net benefit: {n_helped - n_hurt} cases")

print(f"\nResults saved to: {OUTPUT_DIR}")
print("=" * 80)
