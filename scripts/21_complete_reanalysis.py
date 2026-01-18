#!/usr/bin/env python3
"""
Complete Re-analysis with Correct CA19-9 Data for Both Patients and Controls
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
                             recall_score, precision_score, confusion_matrix, f1_score)
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "complete_reanalysis"
OUTPUT_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300

print("=" * 80)
print("COMPLETE PANCREATIC CANCER ANALYSIS - WITH CORRECT CA19-9 DATA")
print("=" * 80)

# ============================================================================
# STEP 1: Load and Integrate All Data
# ============================================================================
print("\n[STEP 1] Loading and integrating all data...")
print("-" * 80)

# Load clinical data
patients_clinical = pd.read_excel(DATA_DIR / "환자군과 대조군260117.xlsx", sheet_name='환자군')
controls_clinical = pd.read_excel(DATA_DIR / "환자군과 대조군260117.xlsx", sheet_name='대조군')

print(f"Clinical data loaded:")
print(f"  Patients: {len(patients_clinical)} (CA19-9: {(~patients_clinical['CA 19-9'].isna()).sum()})")
print(f"  Controls: {len(controls_clinical)} (CA19-9: {(~controls_clinical['CA 19-9'].isna()).sum()})")

# Load methylation data
dmrseq_results = pd.read_csv(RESULTS_DIR / "dss_results/dmrseq/20260117_cut0.05_minR3_bp1000/dmrseq_Significant_DMRs.tsv", sep='\t')
print(f"\nDMRseq results: {len(dmrseq_results)} significant DMRs")

# Load sample methylation data
sample_mapping = pd.read_csv(RESULTS_DIR / "clinical_integration/sample_file_mapping.csv")
print(f"Sample mapping: {len(sample_mapping)} samples")

# Extract methylation levels for each DMR from cov files
print("\nExtracting methylation levels from DMRs...")

def extract_methylation_from_cov(cov_file, dmr_regions):
    """Extract methylation levels for DMR regions"""
    import gzip
    
    methylation_values = []
    
    for idx, dmr in dmr_regions.iterrows():
        chr_name = dmr['seqnames']
        start = dmr['start']
        end = dmr['end']
        
        region_meth = []
        
        try:
            with gzip.open(cov_file, 'rt') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue
                    
                    chrom = parts[0]
                    pos = int(parts[1])
                    
                    if chrom == chr_name and start <= pos <= end:
                        meth_count = int(parts[4])
                        unmeth_count = int(parts[5])
                        total = meth_count + unmeth_count
                        
                        if total > 0:
                            meth_level = (meth_count / total) * 100
                            region_meth.append(meth_level)
        except:
            pass
        
        if len(region_meth) > 0:
            methylation_values.append(np.mean(region_meth))
        else:
            methylation_values.append(np.nan)
    
    return methylation_values

# Process all samples
print("Processing samples...")
all_samples = []

for idx, row in sample_mapping.iterrows():
    sample_id = row['sample_id']
    group = row['group']
    file_path = row['file_path']
    
    if (idx + 1) % 20 == 0:
        print(f"  Processed: {idx + 1}/{len(sample_mapping)}")
    
    meth_values = extract_methylation_from_cov(file_path, dmrseq_results)
    
    sample_data = {
        'sample_id': sample_id,
        'group': group,
        'label': 1 if group == 'patient' else 0
    }
    
    for dmr_idx, meth_val in enumerate(meth_values):
        sample_data[f'DMR_{dmr_idx+1}'] = meth_val
    
    all_samples.append(sample_data)

methylation_df = pd.DataFrame(all_samples)
print(f"\nMethylation data extracted: {methylation_df.shape}")

# Fill missing methylation values with median
dmr_cols_temp = [col for col in methylation_df.columns if col.startswith('DMR_')]
for col in dmr_cols_temp:
    methylation_df[col] = methylation_df[col].fillna(methylation_df[col].median())
print(f"Missing methylation values filled")

# Merge with clinical data
print("\nMerging with clinical data...")

# Prepare clinical data
patients_clinical_clean = patients_clinical[['환자 번호', '나이', '성별', '병기', 'CA 19-9']].copy()
patients_clinical_clean.columns = ['sample_id', 'age', 'sex', 'stage', 'ca19_9']

controls_clinical_clean = controls_clinical[['환자번호', '나이', '성별', 'CA 19-9']].copy()
controls_clinical_clean.columns = ['sample_id', 'age', 'sex', 'ca19_9']
controls_clinical_clean['stage'] = np.nan

# Merge
integrated_df = methylation_df.merge(
    pd.concat([patients_clinical_clean, controls_clinical_clean], ignore_index=True),
    on='sample_id',
    how='left'
)

print(f"Integrated data: {integrated_df.shape}")
print(f"  Patients: {sum(integrated_df['label'] == 1)}")
print(f"  Controls: {sum(integrated_df['label'] == 0)}")
print(f"  CA19-9 missing: {integrated_df['ca19_9'].isna().sum()}")

# Save integrated data
integrated_df.to_csv(OUTPUT_DIR / "integrated_data_complete.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR / 'integrated_data_complete.csv'}")

# ============================================================================
# STEP 2: Overall Model Performance
# ============================================================================
print("\n[STEP 2] Overall model performance...")
print("=" * 80)

dmr_cols = [col for col in integrated_df.columns if col.startswith('DMR_')]
clinical_cols = ['age', 'sex', 'ca19_9']

X_full = integrated_df[dmr_cols + clinical_cols].copy()
y = integrated_df['label']

# Handle missing values
for col in clinical_cols:
    if col == 'sex':
        X_full[col] = X_full[col].fillna(X_full[col].mode()[0] if len(X_full[col].mode()) > 0 else 1)
    else:
        X_full[col] = X_full[col].fillna(X_full[col].median())

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models_results = {}

# Model 1: Methylation only
print("\nModel 1: Methylation only")
X_meth = X_full[dmr_cols].copy()
scaler_meth = StandardScaler()
X_meth_scaled = scaler_meth.fit_transform(X_meth)

lr_meth = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_meth = cross_val_predict(lr_meth, X_meth_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_meth = cross_val_predict(lr_meth, X_meth_scaled, y, cv=cv)

models_results['Methylation only'] = {
    'accuracy': accuracy_score(y, y_pred_meth),
    'auc': roc_auc_score(y, y_prob_meth),
    'sensitivity': recall_score(y, y_pred_meth),
    'specificity': recall_score(y, y_pred_meth, pos_label=0),
    'precision': precision_score(y, y_pred_meth),
    'f1': f1_score(y, y_pred_meth),
    'y_prob': y_prob_meth
}

print(f"  AUC: {models_results['Methylation only']['auc']:.3f}")
print(f"  Accuracy: {models_results['Methylation only']['accuracy']:.3f}")
print(f"  Sensitivity: {models_results['Methylation only']['sensitivity']:.3f}")

# Model 2: Methylation + Age + Sex
print("\nModel 2: Methylation + Age + Sex")
X_no_ca199 = X_full[dmr_cols + ['age', 'sex']].copy()
scaler_no_ca199 = StandardScaler()
X_no_ca199_scaled = scaler_no_ca199.fit_transform(X_no_ca199)

lr_no_ca199 = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_no_ca199 = cross_val_predict(lr_no_ca199, X_no_ca199_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_no_ca199 = cross_val_predict(lr_no_ca199, X_no_ca199_scaled, y, cv=cv)

models_results['Methylation + Age + Sex'] = {
    'accuracy': accuracy_score(y, y_pred_no_ca199),
    'auc': roc_auc_score(y, y_prob_no_ca199),
    'sensitivity': recall_score(y, y_pred_no_ca199),
    'specificity': recall_score(y, y_pred_no_ca199, pos_label=0),
    'precision': precision_score(y, y_pred_no_ca199),
    'f1': f1_score(y, y_pred_no_ca199),
    'y_prob': y_prob_no_ca199
}

print(f"  AUC: {models_results['Methylation + Age + Sex']['auc']:.3f}")
print(f"  Accuracy: {models_results['Methylation + Age + Sex']['accuracy']:.3f}")
print(f"  Sensitivity: {models_results['Methylation + Age + Sex']['sensitivity']:.3f}")

# Model 3: Full model with CA19-9
print("\nModel 3: Methylation + Age + Sex + CA19-9")
X_with_ca199 = X_full[dmr_cols + clinical_cols].copy()
scaler_with_ca199 = StandardScaler()
X_with_ca199_scaled = scaler_with_ca199.fit_transform(X_with_ca199)

lr_with_ca199 = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_with_ca199 = cross_val_predict(lr_with_ca199, X_with_ca199_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_with_ca199 = cross_val_predict(lr_with_ca199, X_with_ca199_scaled, y, cv=cv)

models_results['Methylation + Age + Sex + CA19-9'] = {
    'accuracy': accuracy_score(y, y_pred_with_ca199),
    'auc': roc_auc_score(y, y_prob_with_ca199),
    'sensitivity': recall_score(y, y_pred_with_ca199),
    'specificity': recall_score(y, y_pred_with_ca199, pos_label=0),
    'precision': precision_score(y, y_pred_with_ca199),
    'f1': f1_score(y, y_pred_with_ca199),
    'y_prob': y_prob_with_ca199
}

print(f"  AUC: {models_results['Methylation + Age + Sex + CA19-9']['auc']:.3f}")
print(f"  Accuracy: {models_results['Methylation + Age + Sex + CA19-9']['accuracy']:.3f}")
print(f"  Sensitivity: {models_results['Methylation + Age + Sex + CA19-9']['sensitivity']:.3f}")

# ============================================================================
# STEP 3: Stage-Specific Analysis
# ============================================================================
print("\n[STEP 3] Stage-specific analysis...")
print("=" * 80)

patients_df = integrated_df[integrated_df['label'] == 1].copy()
controls_df = integrated_df[integrated_df['label'] == 0].copy()

early_stage = patients_df[patients_df['stage'].isin([1, 2])].copy()
late_stage = patients_df[patients_df['stage'].isin([3, 4])].copy()

print(f"\nEarly-stage (1-2): {len(early_stage)} patients")
print(f"Late-stage (3-4): {len(late_stage)} patients")
print(f"Controls: {len(controls_df)}")

# Early-stage analysis
print("\nEarly-stage vs Controls:")
early_vs_control = pd.concat([early_stage, controls_df], ignore_index=True)
X_early = early_vs_control[dmr_cols + clinical_cols].copy()
y_early = early_vs_control['label']

for col in clinical_cols:
    if col == 'sex':
        X_early[col] = X_early[col].fillna(X_early[col].mode()[0] if len(X_early[col].mode()) > 0 else 1)
    else:
        X_early[col] = X_early[col].fillna(X_early[col].median())

scaler_early = StandardScaler()
X_early_scaled = scaler_early.fit_transform(X_early)

lr_early = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_early = cross_val_predict(lr_early, X_early_scaled, y_early, cv=cv, method='predict_proba')[:, 1]
y_pred_early = cross_val_predict(lr_early, X_early_scaled, y_early, cv=cv)

acc_early = accuracy_score(y_early, y_pred_early)
auc_early = roc_auc_score(y_early, y_prob_early)
sens_early = recall_score(y_early, y_pred_early)
spec_early = recall_score(y_early, y_pred_early, pos_label=0)

print(f"  AUC: {auc_early:.3f}, Sensitivity: {sens_early:.3f} ({sens_early*100:.1f}%)")
print(f"  Detected: {int(sens_early * len(early_stage))}/{len(early_stage)} early-stage patients")

# Late-stage analysis
print("\nLate-stage vs Controls:")
late_vs_control = pd.concat([late_stage, controls_df], ignore_index=True)
X_late = late_vs_control[dmr_cols + clinical_cols].copy()
y_late = late_vs_control['label']

for col in clinical_cols:
    if col == 'sex':
        X_late[col] = X_late[col].fillna(X_late[col].mode()[0] if len(X_late[col].mode()) > 0 else 1)
    else:
        X_late[col] = X_late[col].fillna(X_late[col].median())

scaler_late = StandardScaler()
X_late_scaled = scaler_late.fit_transform(X_late)

lr_late = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_late = cross_val_predict(lr_late, X_late_scaled, y_late, cv=cv, method='predict_proba')[:, 1]
y_pred_late = cross_val_predict(lr_late, X_late_scaled, y_late, cv=cv)

acc_late = accuracy_score(y_late, y_pred_late)
auc_late = roc_auc_score(y_late, y_prob_late)
sens_late = recall_score(y_late, y_pred_late)
spec_late = recall_score(y_late, y_pred_late, pos_label=0)

print(f"  AUC: {auc_late:.3f}, Sensitivity: {sens_late:.3f} ({sens_late*100:.1f}%)")
print(f"  Detected: {int(sens_late * len(late_stage))}/{len(late_stage)} late-stage patients")

# ============================================================================
# STEP 4: CA19-9 Analysis
# ============================================================================
print("\n[STEP 4] CA19-9 contribution analysis...")
print("=" * 80)

cancer_ca199 = integrated_df[integrated_df['label']==1]['ca19_9'].dropna()
control_ca199 = integrated_df[integrated_df['label']==0]['ca19_9'].dropna()

print(f"\nCA19-9 Statistics:")
print(f"  Cancer: mean={cancer_ca199.mean():.1f}, median={cancer_ca199.median():.1f}")
print(f"  Control: mean={control_ca199.mean():.1f}, median={control_ca199.median():.1f}")
print(f"  Ratio: {cancer_ca199.mean() / control_ca199.mean():.1f}x")

stat, pval = stats.mannwhitneyu(cancer_ca199, control_ca199, alternative='two-sided')
print(f"  Mann-Whitney U test: p = {pval:.2e}")

auc_improvement = models_results['Methylation + Age + Sex + CA19-9']['auc'] - models_results['Methylation + Age + Sex']['auc']
print(f"\nCA19-9 contribution:")
print(f"  AUC improvement: +{auc_improvement:.3f}")

# ============================================================================
# STEP 5: Generate Summary Report
# ============================================================================
print("\n[STEP 5] Generating summary report...")
print("=" * 80)

summary_df = pd.DataFrame({
    'Model': list(models_results.keys()),
    'Accuracy': [m['accuracy'] for m in models_results.values()],
    'AUC': [m['auc'] for m in models_results.values()],
    'Sensitivity': [m['sensitivity'] for m in models_results.values()],
    'Specificity': [m['specificity'] for m in models_results.values()],
    'Precision': [m['precision'] for m in models_results.values()],
    'F1-Score': [m['f1'] for m in models_results.values()]
})

summary_df.to_csv(OUTPUT_DIR / "overall_performance_summary.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR / 'overall_performance_summary.csv'}")

stage_summary_df = pd.DataFrame({
    'Stage': ['Early (1-2)', 'Late (3-4)'],
    'N_Patients': [len(early_stage), len(late_stage)],
    'Accuracy': [acc_early, acc_late],
    'AUC': [auc_early, auc_late],
    'Sensitivity': [sens_early, sens_late],
    'Specificity': [spec_early, spec_late]
})

stage_summary_df.to_csv(OUTPUT_DIR / "stage_specific_summary.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'stage_specific_summary.csv'}")

# ============================================================================
# STEP 6: Visualizations
# ============================================================================
print("\n[STEP 6] Generating visualizations...")

# ROC Curves
fig, ax = plt.subplots(figsize=(10, 8))

for model_name, results in models_results.items():
    fpr, tpr, _ = roc_curve(y, results['y_prob'])
    ax.plot(fpr, tpr, label=f"{model_name} (AUC={results['auc']:.3f})", linewidth=2.5)

ax.plot([0, 1], [0, 1], 'k--', label='Random', linewidth=2)
ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
ax.set_title('ROC Curves: Overall Performance', fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=11, loc='lower right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curves_overall.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'roc_curves_overall.png'}")
plt.close()

# Stage-specific ROC
fig, ax = plt.subplots(figsize=(10, 8))

fpr_early, tpr_early, _ = roc_curve(y_early, y_prob_early)
fpr_late, tpr_late, _ = roc_curve(y_late, y_prob_late)

ax.plot(fpr_early, tpr_early, label=f'Early-stage (1-2) (AUC={auc_early:.3f})', linewidth=3, color='#2ecc71')
ax.plot(fpr_late, tpr_late, label=f'Late-stage (3-4) (AUC={auc_late:.3f})', linewidth=3, color='#e74c3c')
ax.plot([0, 1], [0, 1], 'k--', label='Random', linewidth=2)

ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
ax.set_title('ROC Curves: Stage-Specific Performance', fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=12, loc='lower right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curves_stage_specific.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'roc_curves_stage_specific.png'}")
plt.close()

print("\n" + "=" * 80)
print("COMPLETE REANALYSIS FINISHED!")
print("=" * 80)
print(f"\nResults saved to: {OUTPUT_DIR}")
print("\nSummary:")
print(summary_df.to_string(index=False))
print("\nStage-specific:")
print(stage_summary_df.to_string(index=False))
print("=" * 80)
