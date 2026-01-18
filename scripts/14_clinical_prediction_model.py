#!/usr/bin/env python3
"""
DMRseq 결과와 임상 데이터를 통합하여 췌장암 예측 모델 구축
"""

import pandas as pd
import numpy as np
import gzip
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold, LeaveOneOut
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import (roc_auc_score, roc_curve, classification_report, 
                             confusion_matrix, accuracy_score, precision_recall_curve,
                             precision_score, recall_score, f1_score)
from sklearn.impute import SimpleImputer
import warnings
warnings.filterwarnings('ignore')

# 경로 설정
BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
DMRSEQ_DIR = RESULTS_DIR / "dss_results" / "dmrseq"
OUTPUT_DIR = RESULTS_DIR / "clinical_prediction"
OUTPUT_DIR.mkdir(exist_ok=True)

# 시각화 설정
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

print("=" * 80)
print("췌장암 예측 모델: DMRseq 메틸화 + 임상 데이터 통합 분석")
print("=" * 80)

# 1. 임상 데이터 로드
print("\n[1] 임상 데이터 로드 중...")
clinical_file = DATA_DIR / "환자군과 대조군260117.xlsx"
clinical_df = pd.read_excel(clinical_file)

print(f"임상 데이터: {clinical_df.shape}")
print(f"컬럼: {clinical_df.columns.tolist()}")

# 컬럼명 정리
clinical_df.columns = ['patient_id', 'name', 'age', 'sex', 'stage', 'ca19_9']
clinical_df['group'] = 'patient'

print("\n임상 데이터 통계:")
print(clinical_df[['age', 'sex', 'stage', 'ca19_9']].describe())

# 2. DMRseq 결과 로드
print("\n[2] DMRseq 결과 로드 중...")

# 가장 최근 DMRseq 결과 폴더 찾기
dmrseq_folders = sorted([d for d in DMRSEQ_DIR.glob("*") if d.is_dir()])
if not dmrseq_folders:
    print("ERROR: DMRseq 결과를 찾을 수 없습니다.")
    print(f"경로 확인: {DMRSEQ_DIR}")
    exit(1)

latest_dmrseq = dmrseq_folders[-1]
print(f"DMRseq 결과 폴더: {latest_dmrseq}")

# Significant DMRs 로드
sig_dmr_file = latest_dmrseq / "dmrseq_Significant_DMRs.tsv"
if not sig_dmr_file.exists():
    sig_dmr_file = latest_dmrseq / "dmrseq_DMRs.tsv"

if not sig_dmr_file.exists():
    print("ERROR: DMR 결과 파일을 찾을 수 없습니다.")
    exit(1)

dmrs = pd.read_csv(sig_dmr_file, sep='\t')
print(f"DMRs 로드: {len(dmrs)} 개")
print(f"DMR 컬럼: {dmrs.columns.tolist()}")

# 상위 DMR 선택 (qval 기준)
if 'qval' in dmrs.columns:
    top_dmrs = dmrs.nsmallest(50, 'qval')
else:
    top_dmrs = dmrs.head(50)

print(f"\n상위 {len(top_dmrs)}개 DMR 선택")

# 3. 각 샘플의 메틸화 수준 추출
print("\n[3] 샘플별 메틸화 수준 추출 중...")

patients_dir = DATA_DIR / "patients"
normals_dir = DATA_DIR / "normals"

def extract_patient_id(filename):
    """파일명에서 환자 번호 추출"""
    import re
    basename = filename.name
    if basename.startswith("P-"):
        match = re.match(r"P-(\d+)", basename)
        if match:
            return int(match.group(1))
    match = re.match(r"(\d+)", basename)
    if match:
        return int(match.group(1))
    return None

def extract_normal_id(filename):
    """대조군 파일명에서 번호 추출"""
    import re
    basename = filename.name
    match = re.match(r"(\d+)N", basename)
    if match:
        return int(match.group(1))
    return None

def read_methylation_in_regions(cov_file, regions):
    """특정 영역의 메틸화 수준 계산"""
    try:
        with gzip.open(cov_file, 'rt') as f:
            cov_data = pd.read_csv(f, sep='\t', header=None,
                                   names=['chr', 'start', 'end', 'meth_pct', 'count_m', 'count_u'])
        
        meth_values = []
        for _, region in regions.iterrows():
            chr_name = str(region['seqnames'])
            start = region['start']
            end = region['end']
            
            # 해당 영역의 CpG 사이트 필터링
            region_data = cov_data[
                (cov_data['chr'] == chr_name) &
                (cov_data['start'] >= start) &
                (cov_data['start'] <= end)
            ]
            
            if len(region_data) > 0:
                # 평균 메틸화 비율 계산
                total_m = region_data['count_m'].sum()
                total_u = region_data['count_u'].sum()
                if (total_m + total_u) > 0:
                    meth_pct = total_m / (total_m + total_u) * 100
                else:
                    meth_pct = 0
            else:
                meth_pct = 0
            
            meth_values.append(meth_pct)
        
        return meth_values
    except Exception as e:
        print(f"  오류 ({cov_file.name}): {e}")
        return [0] * len(regions)

# 환자 메틸화 데이터 추출
print("환자 샘플 처리 중...")
patient_meth_data = []
patient_files = list(patients_dir.glob("*.cov.gz"))

for pf in patient_files:
    pid = extract_patient_id(pf)
    if pid is None:
        continue
    
    meth_values = read_methylation_in_regions(pf, top_dmrs)
    
    row_data = {'sample_id': pid, 'group': 'patient'}
    for i, val in enumerate(meth_values):
        row_data[f'DMR_{i+1}'] = val
    
    patient_meth_data.append(row_data)
    
    if len(patient_meth_data) % 10 == 0:
        print(f"  처리 완료: {len(patient_meth_data)}/{len(patient_files)}")

# 대조군 메틸화 데이터 추출
print("대조군 샘플 처리 중...")
normal_meth_data = []
normal_files = list(normals_dir.glob("*.cov.gz"))

for nf in normal_files:
    nid = extract_normal_id(nf)
    if nid is None:
        continue
    
    meth_values = read_methylation_in_regions(nf, top_dmrs)
    
    row_data = {'sample_id': nid, 'group': 'normal'}
    for i, val in enumerate(meth_values):
        row_data[f'DMR_{i+1}'] = val
    
    normal_meth_data.append(row_data)
    
    if len(normal_meth_data) % 10 == 0:
        print(f"  처리 완료: {len(normal_meth_data)}/{len(normal_files)}")

# 메틸화 데이터프레임 생성
meth_df = pd.DataFrame(patient_meth_data + normal_meth_data)
print(f"\n메틸화 데이터: {meth_df.shape}")

# 4. 임상 데이터와 메틸화 데이터 통합
print("\n[4] 데이터 통합 중...")

# 환자군만 임상 데이터와 매칭
patient_meth = meth_df[meth_df['group'] == 'patient'].copy()
patient_integrated = patient_meth.merge(
    clinical_df[['patient_id', 'age', 'sex', 'stage', 'ca19_9']], 
    left_on='sample_id', 
    right_on='patient_id', 
    how='left'
)
patient_integrated = patient_integrated.drop('patient_id', axis=1)

# 대조군 (임상 데이터 없음)
normal_meth = meth_df[meth_df['group'] == 'normal'].copy()
normal_meth['age'] = np.nan
normal_meth['sex'] = np.nan
normal_meth['stage'] = np.nan
normal_meth['ca19_9'] = np.nan

# 전체 데이터 합치기
full_data = pd.concat([patient_integrated, normal_meth], ignore_index=True)
full_data['label'] = (full_data['group'] == 'patient').astype(int)

print(f"통합 데이터: {full_data.shape}")
print(f"  환자: {sum(full_data['label'] == 1)}")
print(f"  대조군: {sum(full_data['label'] == 0)}")

# 데이터 저장
full_data.to_csv(OUTPUT_DIR / "integrated_data.csv", index=False)
print(f"\n통합 데이터 저장: {OUTPUT_DIR / 'integrated_data.csv'}")

# 5. 예측 모델 구축
print("\n[5] 예측 모델 구축 중...")

# 특징 선택
dmr_cols = [col for col in full_data.columns if col.startswith('DMR_')]
clinical_cols = ['age', 'sex', 'stage', 'ca19_9']

# 모델 1: 메틸화 데이터만
print("\n=== 모델 1: 메틸화 특징만 ===")
X_meth = full_data[dmr_cols].fillna(0)
y = full_data['label']

X_train, X_test, y_train, y_test = train_test_split(X_meth, y, test_size=0.3, random_state=42, stratify=y)

scaler_meth = StandardScaler()
X_train_scaled = scaler_meth.fit_transform(X_train)
X_test_scaled = scaler_meth.transform(X_test)

# Logistic Regression
lr_meth = LogisticRegression(max_iter=1000, random_state=42)
lr_meth.fit(X_train_scaled, y_train)
y_pred_meth = lr_meth.predict(X_test_scaled)
y_prob_meth = lr_meth.predict_proba(X_test_scaled)[:, 1]

print(f"정확도: {accuracy_score(y_test, y_pred_meth):.3f}")
print(f"AUC: {roc_auc_score(y_test, y_prob_meth):.3f}")
print(f"민감도: {recall_score(y_test, y_pred_meth):.3f}")
print(f"특이도: {recall_score(y_test, y_pred_meth, pos_label=0):.3f}")

# 모델 2: 메틸화 + 임상 데이터 (환자군만)
print("\n=== 모델 2: 메틸화 + 임상 특징 (환자군만) ===")

# 환자군 데이터만 사용
patient_data = full_data[full_data['label'] == 1].copy()

# 결측치 처리
for col in clinical_cols:
    if col == 'sex':
        patient_data[col] = patient_data[col].fillna(patient_data[col].mode()[0] if len(patient_data[col].mode()) > 0 else 1)
    else:
        patient_data[col] = patient_data[col].fillna(patient_data[col].median())

# 대조군 데이터 (임상 특징은 평균값 사용)
normal_data = full_data[full_data['label'] == 0].copy()
for col in clinical_cols:
    if col == 'sex':
        normal_data[col] = patient_data[col].mode()[0] if len(patient_data[col].mode()) > 0 else 1
    else:
        normal_data[col] = patient_data[col].median()

# 데이터 합치기
combined_data = pd.concat([patient_data, normal_data], ignore_index=True)

X_combined = combined_data[dmr_cols + clinical_cols]
y_combined = combined_data['label']

X_train_c, X_test_c, y_train_c, y_test_c = train_test_split(
    X_combined, y_combined, test_size=0.3, random_state=42, stratify=y_combined
)

scaler_combined = StandardScaler()
X_train_c_scaled = scaler_combined.fit_transform(X_train_c)
X_test_c_scaled = scaler_combined.transform(X_test_c)

# Logistic Regression
lr_combined = LogisticRegression(max_iter=1000, random_state=42)
lr_combined.fit(X_train_c_scaled, y_train_c)
y_pred_combined = lr_combined.predict(X_test_c_scaled)
y_prob_combined = lr_combined.predict_proba(X_test_c_scaled)[:, 1]

print(f"정확도: {accuracy_score(y_test_c, y_pred_combined):.3f}")
print(f"AUC: {roc_auc_score(y_test_c, y_prob_combined):.3f}")
print(f"민감도: {recall_score(y_test_c, y_pred_combined):.3f}")
print(f"특이도: {recall_score(y_test_c, y_pred_combined, pos_label=0):.3f}")

# 6. 시각화
print("\n[6] 시각화 생성 중...")

# ROC Curve 비교
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 모델 1: 메틸화만
fpr_meth, tpr_meth, _ = roc_curve(y_test, y_prob_meth)
auc_meth = roc_auc_score(y_test, y_prob_meth)
axes[0].plot(fpr_meth, tpr_meth, label=f'메틸화만 (AUC={auc_meth:.3f})', linewidth=2)
axes[0].plot([0, 1], [0, 1], 'k--', label='Random')
axes[0].set_xlabel('False Positive Rate', fontsize=12)
axes[0].set_ylabel('True Positive Rate', fontsize=12)
axes[0].set_title('ROC Curve - 메틸화 특징만', fontsize=14, fontweight='bold')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# 모델 2: 메틸화 + 임상
fpr_comb, tpr_comb, _ = roc_curve(y_test_c, y_prob_combined)
auc_comb = roc_auc_score(y_test_c, y_prob_combined)
axes[1].plot(fpr_comb, tpr_comb, label=f'메틸화+임상 (AUC={auc_comb:.3f})', linewidth=2, color='red')
axes[1].plot([0, 1], [0, 1], 'k--', label='Random')
axes[1].set_xlabel('False Positive Rate', fontsize=12)
axes[1].set_ylabel('True Positive Rate', fontsize=12)
axes[1].set_title('ROC Curve - 메틸화 + 임상 특징', fontsize=14, fontweight='bold')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "roc_curves_comparison.png", dpi=300, bbox_inches='tight')
print(f"저장: {OUTPUT_DIR / 'roc_curves_comparison.png'}")

# Feature Importance (메틸화+임상 모델)
fig, ax = plt.subplots(figsize=(10, 8))
feature_names = dmr_cols + clinical_cols
importances = np.abs(lr_combined.coef_[0])
indices = np.argsort(importances)[-20:]  # 상위 20개

ax.barh(range(len(indices)), importances[indices])
ax.set_yticks(range(len(indices)))
ax.set_yticklabels([feature_names[i] for i in indices])
ax.set_xlabel('Absolute Coefficient', fontsize=12)
ax.set_title('상위 20개 중요 특징 (메틸화 + 임상)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "feature_importance.png", dpi=300, bbox_inches='tight')
print(f"저장: {OUTPUT_DIR / 'feature_importance.png'}")

# 결과 요약
summary = {
    '모델': ['메틸화만', '메틸화+임상'],
    '정확도': [
        accuracy_score(y_test, y_pred_meth),
        accuracy_score(y_test_c, y_pred_combined)
    ],
    'AUC': [
        roc_auc_score(y_test, y_prob_meth),
        roc_auc_score(y_test_c, y_prob_combined)
    ],
    '민감도': [
        recall_score(y_test, y_pred_meth),
        recall_score(y_test_c, y_pred_combined)
    ],
    '특이도': [
        recall_score(y_test, y_pred_meth, pos_label=0),
        recall_score(y_test_c, y_pred_combined, pos_label=0)
    ]
}

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUTPUT_DIR / "model_performance_summary.csv", index=False)
print(f"\n저장: {OUTPUT_DIR / 'model_performance_summary.csv'}")

print("\n" + "=" * 80)
print("분석 완료!")
print("=" * 80)
print(f"\n결과 저장 위치: {OUTPUT_DIR}")
print("\n모델 성능 요약:")
print(summary_df.to_string(index=False))
