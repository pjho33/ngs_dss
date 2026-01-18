#!/usr/bin/env python3
"""
병기별 췌장암 진단 성능 분석 - 특히 초기 암(1-2기) 진단율 평가
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

# 경로 설정
BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "stage_specific_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

# 시각화 설정
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['figure.dpi'] = 300

print("=" * 80)
print("병기별 췌장암 진단 성능 분석")
print("=" * 80)

# 1. 통합 데이터 로드
print("\n[1] 통합 데이터 로드 중...")
data_file = RESULTS_DIR / "clinical_prediction" / "integrated_data.csv"
df = pd.read_csv(data_file)

print(f"전체 데이터: {df.shape}")
print(f"  환자: {sum(df['label'] == 1)}")
print(f"  대조군: {sum(df['label'] == 0)}")

# 2. 병기별 데이터 분석
print("\n[2] 병기별 데이터 분석...")

# 환자군만 추출
patient_df = df[df['label'] == 1].copy()

print(f"\n환자군 병기 분포:")
stage_counts = patient_df['stage'].value_counts().sort_index()
print(stage_counts)

# 결측치 확인
print(f"\n병기 결측치: {patient_df['stage'].isna().sum()}명")

# 병기별 그룹 생성
patient_df['stage_group'] = patient_df['stage'].apply(
    lambda x: 'Early (1-2기)' if x in [1, 2] else ('Late (3-4기)' if x in [3, 4] else 'Unknown')
)

print(f"\n병기 그룹 분포:")
print(patient_df['stage_group'].value_counts())

# 3. 초기 암(1-2기) vs 대조군 분석
print("\n[3] 초기 암(1-2기) vs 대조군 진단 성능 분석...")

# 초기 암 환자
early_stage = patient_df[patient_df['stage'].isin([1, 2])].copy()
print(f"\n초기 암 환자: {len(early_stage)}명")

# 후기 암 환자
late_stage = patient_df[patient_df['stage'].isin([3, 4])].copy()
print(f"후기 암 환자: {len(late_stage)}명")

# 대조군
normal_df = df[df['label'] == 0].copy()
print(f"대조군: {len(normal_df)}명")

# DMR 특징 컬럼
dmr_cols = [col for col in df.columns if col.startswith('DMR_')]
clinical_cols = ['age', 'sex', 'ca19_9']  # stage는 제외 (예측 대상이므로)

print(f"\nDMR 특징 수: {len(dmr_cols)}")

# 4. 모델 1: 초기 암 vs 대조군
print("\n" + "=" * 80)
print("모델 1: 초기 암(1-2기) vs 대조군")
print("=" * 80)

early_vs_normal = pd.concat([early_stage, normal_df], ignore_index=True)
X_early = early_vs_normal[dmr_cols + clinical_cols].fillna(0)
y_early = early_vs_normal['label']

# 결측치 처리
for col in clinical_cols:
    if col in ['sex']:
        X_early[col] = X_early[col].fillna(X_early[col].mode()[0] if len(X_early[col].mode()) > 0 else 1)
    else:
        X_early[col] = X_early[col].fillna(X_early[col].median())

# Cross-validation으로 평가
scaler_early = StandardScaler()
X_early_scaled = scaler_early.fit_transform(X_early)

lr_early = LogisticRegression(max_iter=1000, random_state=42, class_weight='balanced')
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Cross-validation 예측
y_pred_early_cv = cross_val_predict(lr_early, X_early_scaled, y_early, cv=cv)
y_prob_early_cv = cross_val_predict(lr_early, X_early_scaled, y_early, cv=cv, method='predict_proba')[:, 1]

# 성능 평가
acc_early = accuracy_score(y_early, y_pred_early_cv)
auc_early = roc_auc_score(y_early, y_prob_early_cv)
sens_early = recall_score(y_early, y_pred_early_cv)
spec_early = recall_score(y_early, y_pred_early_cv, pos_label=0)
prec_early = precision_score(y_early, y_pred_early_cv)
f1_early = f1_score(y_early, y_pred_early_cv)

print(f"\n정확도: {acc_early:.3f}")
print(f"AUC: {auc_early:.3f}")
print(f"민감도 (초기 암 탐지율): {sens_early:.3f}")
print(f"특이도: {spec_early:.3f}")
print(f"정밀도: {prec_early:.3f}")
print(f"F1-score: {f1_early:.3f}")

# Confusion Matrix
cm_early = confusion_matrix(y_early, y_pred_early_cv)
print(f"\nConfusion Matrix:")
print(f"  TN: {cm_early[0,0]}, FP: {cm_early[0,1]}")
print(f"  FN: {cm_early[1,0]}, TP: {cm_early[1,1]}")

# 5. 모델 2: 후기 암 vs 대조군
print("\n" + "=" * 80)
print("모델 2: 후기 암(3-4기) vs 대조군")
print("=" * 80)

late_vs_normal = pd.concat([late_stage, normal_df], ignore_index=True)
X_late = late_vs_normal[dmr_cols + clinical_cols].fillna(0)
y_late = late_vs_normal['label']

# 결측치 처리
for col in clinical_cols:
    if col in ['sex']:
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

print(f"\n정확도: {acc_late:.3f}")
print(f"AUC: {auc_late:.3f}")
print(f"민감도 (후기 암 탐지율): {sens_late:.3f}")
print(f"특이도: {spec_late:.3f}")
print(f"정밀도: {prec_late:.3f}")
print(f"F1-score: {f1_late:.3f}")

cm_late = confusion_matrix(y_late, y_pred_late_cv)
print(f"\nConfusion Matrix:")
print(f"  TN: {cm_late[0,0]}, FP: {cm_late[0,1]}")
print(f"  FN: {cm_late[1,0]}, TP: {cm_late[1,1]}")

# 6. 시각화
print("\n[6] 시각화 생성 중...")

# 6-1. ROC Curve 비교
fig, ax = plt.subplots(figsize=(10, 8))

fpr_early, tpr_early, _ = roc_curve(y_early, y_prob_early_cv)
fpr_late, tpr_late, _ = roc_curve(y_late, y_prob_late_cv)

ax.plot(fpr_early, tpr_early, label=f'초기 암(1-2기) vs 대조군 (AUC={auc_early:.3f})', 
        linewidth=3, color='#2ecc71')
ax.plot(fpr_late, tpr_late, label=f'후기 암(3-4기) vs 대조군 (AUC={auc_late:.3f})', 
        linewidth=3, color='#e74c3c')
ax.plot([0, 1], [0, 1], 'k--', label='Random', linewidth=2)

ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
ax.set_title('ROC Curve: 병기별 진단 성능 비교', fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=12, loc='lower right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curve_by_stage.png", dpi=300, bbox_inches='tight')
print(f"저장: {OUTPUT_DIR / 'roc_curve_by_stage.png'}")
plt.close()

# 6-2. 성능 비교 바 차트
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

metrics = ['정확도', 'AUC', '민감도', '특이도']
early_scores = [acc_early, auc_early, sens_early, spec_early]
late_scores = [acc_late, auc_late, sens_late, spec_late]

for idx, (ax, metric, early_val, late_val) in enumerate(zip(axes.flat, metrics, early_scores, late_scores)):
    x = ['초기 암\n(1-2기)', '후기 암\n(3-4기)']
    values = [early_val, late_val]
    colors = ['#2ecc71', '#e74c3c']
    
    bars = ax.bar(x, values, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
    ax.set_ylabel(metric, fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1.1])
    ax.grid(True, alpha=0.3, axis='y')
    
    # 값 표시
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{val:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.suptitle('병기별 진단 성능 비교', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "performance_by_stage.png", dpi=300, bbox_inches='tight')
print(f"저장: {OUTPUT_DIR / 'performance_by_stage.png'}")
plt.close()

# 6-3. Confusion Matrix 시각화
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 초기 암
sns.heatmap(cm_early, annot=True, fmt='d', cmap='Greens', ax=axes[0], 
            cbar_kws={'label': 'Count'}, annot_kws={'size': 14, 'weight': 'bold'})
axes[0].set_title(f'초기 암(1-2기) vs 대조군\n민감도: {sens_early:.1%}', 
                  fontsize=14, fontweight='bold')
axes[0].set_ylabel('실제', fontsize=12, fontweight='bold')
axes[0].set_xlabel('예측', fontsize=12, fontweight='bold')
axes[0].set_xticklabels(['대조군', '환자'])
axes[0].set_yticklabels(['대조군', '환자'])

# 후기 암
sns.heatmap(cm_late, annot=True, fmt='d', cmap='Reds', ax=axes[1], 
            cbar_kws={'label': 'Count'}, annot_kws={'size': 14, 'weight': 'bold'})
axes[1].set_title(f'후기 암(3-4기) vs 대조군\n민감도: {sens_late:.1%}', 
                  fontsize=14, fontweight='bold')
axes[1].set_ylabel('실제', fontsize=12, fontweight='bold')
axes[1].set_xlabel('예측', fontsize=12, fontweight='bold')
axes[1].set_xticklabels(['대조군', '환자'])
axes[1].set_yticklabels(['대조군', '환자'])

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "confusion_matrix_by_stage.png", dpi=300, bbox_inches='tight')
print(f"저장: {OUTPUT_DIR / 'confusion_matrix_by_stage.png'}")
plt.close()

# 7. 결과 요약 저장
print("\n[7] 결과 요약 저장...")

summary = pd.DataFrame({
    '병기': ['초기 암 (1-2기)', '후기 암 (3-4기)'],
    '환자 수': [len(early_stage), len(late_stage)],
    '정확도': [acc_early, acc_late],
    'AUC': [auc_early, auc_late],
    '민감도': [sens_early, sens_late],
    '특이도': [spec_early, spec_late],
    '정밀도': [prec_early, prec_late],
    'F1-score': [f1_early, f1_late]
})

summary.to_csv(OUTPUT_DIR / "stage_specific_performance.csv", index=False)
print(f"저장: {OUTPUT_DIR / 'stage_specific_performance.csv'}")

# 상세 결과
detailed_results = {
    '초기 암 (1-2기)': {
        '환자 수': len(early_stage),
        '대조군 수': len(normal_df),
        '정확도': acc_early,
        'AUC': auc_early,
        '민감도 (초기 암 탐지율)': sens_early,
        '특이도': spec_early,
        '정밀도': prec_early,
        'F1-score': f1_early,
        'True Positive': int(cm_early[1,1]),
        'False Negative': int(cm_early[1,0]),
        'True Negative': int(cm_early[0,0]),
        'False Positive': int(cm_early[0,1])
    },
    '후기 암 (3-4기)': {
        '환자 수': len(late_stage),
        '대조군 수': len(normal_df),
        '정확도': acc_late,
        'AUC': auc_late,
        '민감도 (후기 암 탐지율)': sens_late,
        '특이도': spec_late,
        '정밀도': prec_late,
        'F1-score': f1_late,
        'True Positive': int(cm_late[1,1]),
        'False Negative': int(cm_late[1,0]),
        'True Negative': int(cm_late[0,0]),
        'False Positive': int(cm_late[0,1])
    }
}

detailed_df = pd.DataFrame(detailed_results).T
detailed_df.to_csv(OUTPUT_DIR / "detailed_stage_results.csv")
print(f"저장: {OUTPUT_DIR / 'detailed_stage_results.csv'}")

print("\n" + "=" * 80)
print("병기별 분석 완료!")
print("=" * 80)

print("\n📊 주요 결과:")
print(f"\n초기 암(1-2기) 진단 성능:")
print(f"  - 환자 수: {len(early_stage)}명")
print(f"  - 민감도 (초기 암 탐지율): {sens_early:.1%}")
print(f"  - AUC: {auc_early:.3f}")
print(f"  - 정확도: {acc_early:.1%}")

print(f"\n후기 암(3-4기) 진단 성능:")
print(f"  - 환자 수: {len(late_stage)}명")
print(f"  - 민감도 (후기 암 탐지율): {sens_late:.1%}")
print(f"  - AUC: {auc_late:.3f}")
print(f"  - 정확도: {acc_late:.1%}")

print(f"\n결과 저장 위치: {OUTPUT_DIR}")
