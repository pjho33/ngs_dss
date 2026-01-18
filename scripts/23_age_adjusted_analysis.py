#!/usr/bin/env python3
"""
나이 보정 분석 - Age-Adjusted Analysis
나이 차이로 인한 편향을 제거하고 순수한 메틸화 효과 분석
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
                             recall_score, precision_score, f1_score)
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "age_adjusted_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300

print("=" * 80)
print("나이 보정 분석 (Age-Adjusted Analysis)")
print("=" * 80)

# Load data
integrated = pd.read_csv(RESULTS_DIR / "complete_reanalysis" / "integrated_data_complete.csv")

patients = integrated[integrated['label'] == 1].copy()
controls = integrated[integrated['label'] == 0].copy()

print(f"\n원본 데이터:")
print(f"  환자군: {len(patients)}명, 평균 나이 {patients['age'].mean():.1f}세")
print(f"  대조군: {len(controls)}명, 평균 나이 {controls['age'].mean():.1f}세")
print(f"  나이 차이: {patients['age'].mean() - controls['age'].mean():.1f}세")

# ============================================================================
# 방법 1: 나이를 모델에서 제외 (메틸화만 사용)
# ============================================================================
print("\n" + "=" * 80)
print("방법 1: 나이 제외 분석 (메틸화 + 성별 + CA19-9)")
print("=" * 80)

dmr_cols = [col for col in integrated.columns if col.startswith('DMR_')]
X_no_age = integrated[dmr_cols + ['sex', 'ca19_9']].copy()
y = integrated['label']

# Fill missing values
X_no_age['sex'] = X_no_age['sex'].fillna(X_no_age['sex'].mode()[0])
X_no_age['ca19_9'] = X_no_age['ca19_9'].fillna(X_no_age['ca19_9'].median())

scaler = StandardScaler()
X_no_age_scaled = scaler.fit_transform(X_no_age)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
lr = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')

y_prob_no_age = cross_val_predict(lr, X_no_age_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_no_age = cross_val_predict(lr, X_no_age_scaled, y, cv=cv)

acc_no_age = accuracy_score(y, y_pred_no_age)
auc_no_age = roc_auc_score(y, y_prob_no_age)
sens_no_age = recall_score(y, y_pred_no_age)
spec_no_age = recall_score(y, y_pred_no_age, pos_label=0)

print(f"\n성능 (나이 제외):")
print(f"  AUC: {auc_no_age:.3f}")
print(f"  정확도: {acc_no_age:.3f} ({acc_no_age*100:.1f}%)")
print(f"  민감도: {sens_no_age:.3f} ({sens_no_age*100:.1f}%)")
print(f"  특이도: {spec_no_age:.3f} ({spec_no_age*100:.1f}%)")

# ============================================================================
# 방법 2: 메틸화만 사용 (순수 메틸화 효과)
# ============================================================================
print("\n" + "=" * 80)
print("방법 2: 메틸화만 사용 (순수 메틸화 효과)")
print("=" * 80)

X_meth_only = integrated[dmr_cols].copy()
scaler_meth = StandardScaler()
X_meth_only_scaled = scaler_meth.fit_transform(X_meth_only)

lr_meth = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_meth = cross_val_predict(lr_meth, X_meth_only_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_meth = cross_val_predict(lr_meth, X_meth_only_scaled, y, cv=cv)

acc_meth = accuracy_score(y, y_pred_meth)
auc_meth = roc_auc_score(y, y_prob_meth)
sens_meth = recall_score(y, y_pred_meth)
spec_meth = recall_score(y, y_pred_meth, pos_label=0)

print(f"\n성능 (메틸화만):")
print(f"  AUC: {auc_meth:.3f}")
print(f"  정확도: {acc_meth:.3f} ({acc_meth*100:.1f}%)")
print(f"  민감도: {sens_meth:.3f} ({sens_meth*100:.1f}%)")
print(f"  특이도: {spec_meth:.3f} ({spec_meth*100:.1f}%)")

# ============================================================================
# 방법 3: 나이의 기여도 분석
# ============================================================================
print("\n" + "=" * 80)
print("방법 3: 나이의 기여도 분석")
print("=" * 80)

# Age only
X_age_only = integrated[['age']].copy()
X_age_only = X_age_only.fillna(X_age_only.median())
scaler_age = StandardScaler()
X_age_only_scaled = scaler_age.fit_transform(X_age_only)

lr_age = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
y_prob_age = cross_val_predict(lr_age, X_age_only_scaled, y, cv=cv, method='predict_proba')[:, 1]
y_pred_age = cross_val_predict(lr_age, X_age_only_scaled, y, cv=cv)

acc_age = accuracy_score(y, y_pred_age)
auc_age = roc_auc_score(y, y_prob_age)
sens_age = recall_score(y, y_pred_age)
spec_age = recall_score(y, y_pred_age, pos_label=0)

print(f"\n성능 (나이만):")
print(f"  AUC: {auc_age:.3f}")
print(f"  정확도: {acc_age:.3f} ({acc_age*100:.1f}%)")
print(f"  민감도: {sens_age:.3f} ({sens_age*100:.1f}%)")
print(f"  특이도: {spec_age:.3f} ({spec_age*100:.1f}%)")

print(f"\n⚠️  나이만으로도 AUC {auc_age:.3f} 달성!")
print(f"   이는 나이 차이가 모델 성능에 큰 영향을 미침을 의미합니다.")

# ============================================================================
# 방법 4: 나이 매칭된 서브셋 분석
# ============================================================================
print("\n" + "=" * 80)
print("방법 4: 나이 매칭 분석 (40세 이상 대조군만 사용)")
print("=" * 80)

# 40세 이상 대조군만 선택
controls_40plus = controls[controls['age'] >= 40].copy()
print(f"\n40세 이상 대조군: {len(controls_40plus)}명")
print(f"  평균 나이: {controls_40plus['age'].mean():.1f}세")

if len(controls_40plus) > 0:
    # 나이 매칭된 데이터셋
    matched_data = pd.concat([patients, controls_40plus], ignore_index=True)
    
    X_matched = matched_data[dmr_cols + ['age', 'sex', 'ca19_9']].copy()
    y_matched = matched_data['label']
    
    # Fill missing
    X_matched['sex'] = X_matched['sex'].fillna(X_matched['sex'].mode()[0])
    X_matched['ca19_9'] = X_matched['ca19_9'].fillna(X_matched['ca19_9'].median())
    X_matched['age'] = X_matched['age'].fillna(X_matched['age'].median())
    
    scaler_matched = StandardScaler()
    X_matched_scaled = scaler_matched.fit_transform(X_matched)
    
    lr_matched = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
    y_prob_matched = cross_val_predict(lr_matched, X_matched_scaled, y_matched, cv=cv, method='predict_proba')[:, 1]
    y_pred_matched = cross_val_predict(lr_matched, X_matched_scaled, y_matched, cv=cv)
    
    acc_matched = accuracy_score(y_matched, y_pred_matched)
    auc_matched = roc_auc_score(y_matched, y_prob_matched)
    sens_matched = recall_score(y_matched, y_pred_matched)
    spec_matched = recall_score(y_matched, y_pred_matched, pos_label=0)
    
    print(f"\n성능 (나이 매칭):")
    print(f"  AUC: {auc_matched:.3f}")
    print(f"  정확도: {acc_matched:.3f} ({acc_matched*100:.1f}%)")
    print(f"  민감도: {sens_matched:.3f} ({sens_matched*100:.1f}%)")
    print(f"  특이도: {spec_matched:.3f} ({spec_matched*100:.1f}%)")
    
    print(f"\n나이 차이:")
    print(f"  환자군: {patients['age'].mean():.1f}세")
    print(f"  대조군 (40+): {controls_40plus['age'].mean():.1f}세")
    print(f"  차이: {patients['age'].mean() - controls_40plus['age'].mean():.1f}세")
else:
    print("40세 이상 대조군이 없습니다.")
    auc_matched = None

# ============================================================================
# 요약 및 시각화
# ============================================================================
print("\n" + "=" * 80)
print("분석 결과 요약")
print("=" * 80)

summary_data = {
    '모델': [
        '나이만',
        '메틸화만',
        '메틸화 + 성별 + CA19-9 (나이 제외)',
        '메틸화 + 나이 + 성별 + CA19-9 (원본)',
    ],
    'AUC': [auc_age, auc_meth, auc_no_age, 0.978],
    '정확도': [acc_age, acc_meth, acc_no_age, 0.949],
    '민감도': [sens_age, sens_meth, sens_no_age, 0.942],
    '특이도': [spec_age, spec_meth, spec_no_age, 0.960]
}

if auc_matched is not None:
    summary_data['모델'].append('나이 매칭 (40세 이상)')
    summary_data['AUC'].append(auc_matched)
    summary_data['정확도'].append(acc_matched)
    summary_data['민감도'].append(sens_matched)
    summary_data['특이도'].append(spec_matched)

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(OUTPUT_DIR / "age_adjusted_summary.csv", index=False)

print("\n" + summary_df.to_string(index=False))

# Visualizations
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Age distribution
axes[0].hist(patients['age'], bins=15, alpha=0.6, label=f'환자군 (n={len(patients)})', 
             color='red', edgecolor='black')
axes[0].hist(controls['age'], bins=15, alpha=0.6, label=f'대조군 (n={len(controls)})', 
             color='blue', edgecolor='black')
axes[0].axvline(patients['age'].mean(), color='red', linestyle='--', linewidth=2, 
                label=f'환자 평균: {patients["age"].mean():.1f}세')
axes[0].axvline(controls['age'].mean(), color='blue', linestyle='--', linewidth=2,
                label=f'대조군 평균: {controls["age"].mean():.1f}세')
axes[0].set_xlabel('Age (years)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Frequency', fontsize=12, fontweight='bold')
axes[0].set_title('Age Distribution: Patients vs Controls', fontsize=14, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.3)

# ROC curves comparison
fpr_age, tpr_age, _ = roc_curve(y, y_prob_age)
fpr_meth, tpr_meth, _ = roc_curve(y, y_prob_meth)
fpr_no_age, tpr_no_age, _ = roc_curve(y, y_prob_no_age)

axes[1].plot(fpr_age, tpr_age, label=f'Age only (AUC={auc_age:.3f})', 
             linewidth=2.5, color='#e74c3c', linestyle='--')
axes[1].plot(fpr_meth, tpr_meth, label=f'Methylation only (AUC={auc_meth:.3f})', 
             linewidth=2.5, color='#3498db')
axes[1].plot(fpr_no_age, tpr_no_age, label=f'Methylation + Sex + CA19-9 (AUC={auc_no_age:.3f})', 
             linewidth=3, color='#2ecc71')
axes[1].plot([0, 1], [0, 1], 'k--', linewidth=2)

axes[1].set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
axes[1].set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
axes[1].set_title('ROC Curves: Age Effect Analysis', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=10, loc='lower right')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "age_effect_analysis.png", dpi=300, bbox_inches='tight')
print(f"\nSaved: {OUTPUT_DIR / 'age_effect_analysis.png'}")
plt.close()

# ============================================================================
# 결론
# ============================================================================
print("\n" + "=" * 80)
print("결론 및 권장사항")
print("=" * 80)

print(f"\n1. 나이 효과:")
print(f"   - 나이만으로 AUC {auc_age:.3f} 달성")
print(f"   - 환자와 대조군의 평균 나이 차이: {patients['age'].mean() - controls['age'].mean():.1f}세")
print(f"   - 이는 모델이 암보다는 나이를 학습했을 가능성을 시사합니다.")

print(f"\n2. 순수 메틸화 효과:")
print(f"   - 메틸화만 사용 시 AUC {auc_meth:.3f}")
print(f"   - 나이를 제외하면 AUC {auc_no_age:.3f}")
print(f"   - 메틸화 데이터 자체는 진단 가치가 있습니다.")

print(f"\n3. 권장사항:")
print(f"   ⚠️  현재 대조군의 나이가 너무 젊습니다 (평균 34세 vs 환자 73.5세)")
print(f"   ✓  나이 매칭된 대조군 (60세 이상)을 추가로 확보하는 것이 필요합니다.")
print(f"   ✓  또는 나이를 보정한 분석 결과를 보고해야 합니다.")
print(f"   ✓  현재 결과는 '나이 효과'가 크게 포함되어 있음을 명시해야 합니다.")

print("\n" + "=" * 80)
print(f"Results saved to: {OUTPUT_DIR}")
print("=" * 80)
