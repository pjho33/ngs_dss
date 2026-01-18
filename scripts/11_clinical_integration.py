#!/usr/bin/env python3
"""
EMR-seq 메틸화 데이터와 임상 데이터를 통합하여 췌장암 예측 모델 구축
"""

import pandas as pd
import numpy as np
import os
import re
import gzip
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import (roc_auc_score, roc_curve, classification_report, 
                             confusion_matrix, accuracy_score, precision_recall_curve)
import warnings
warnings.filterwarnings('ignore')

# 경로 설정
BASE_DIR = Path("/home/pjho3/projects/ngs_dss")
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# 출력 디렉토리 생성
OUTPUT_DIR = RESULTS_DIR / "clinical_integration"
OUTPUT_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("EMR-seq 메틸화 데이터 + 임상 데이터 통합 분석")
print("=" * 80)

# 1. 임상 데이터 로드
print("\n[1] 임상 데이터 로드 중...")
clinical_file = DATA_DIR / "환자군과 대조군260117.xlsx"
clinical_df = pd.read_excel(clinical_file)

print(f"임상 데이터 shape: {clinical_df.shape}")
print(f"컬럼: {clinical_df.columns.tolist()}")
print("\n임상 데이터 미리보기:")
print(clinical_df.head(10))

# 2. 환자 및 대조군 파일 매칭
print("\n[2] 환자 및 대조군 파일 매칭 중...")

patients_dir = DATA_DIR / "patients"
normals_dir = DATA_DIR / "normals"

# 환자 파일 리스트
patient_files = list(patients_dir.glob("*.cov.gz"))
normal_files = list(normals_dir.glob("*.cov.gz"))

print(f"환자 파일 수: {len(patient_files)}")
print(f"대조군 파일 수: {len(normal_files)}")

# 환자 번호 추출 함수
def extract_patient_id(filename):
    """파일명에서 환자 번호 추출"""
    basename = os.path.basename(filename)
    
    # P-로 시작하는 경우 (예: P-106-PAN)
    if basename.startswith("P-"):
        match = re.match(r"P-(\d+)", basename)
        if match:
            return int(match.group(1))
    
    # 숫자로 시작하는 경우 (예: 167_S53)
    match = re.match(r"(\d+)", basename)
    if match:
        return int(match.group(1))
    
    return None

def extract_normal_id(filename):
    """대조군 파일명에서 번호 추출 (예: 11N_S10 -> 11)"""
    basename = os.path.basename(filename)
    match = re.match(r"(\d+)N", basename)
    if match:
        return int(match.group(1))
    return None

# 파일 매칭 테이블 생성
file_mapping = []

for pf in patient_files:
    pid = extract_patient_id(str(pf))
    if pid:
        file_mapping.append({
            'sample_id': pid,
            'group': 'patient',
            'file_path': str(pf)
        })

for nf in normal_files:
    nid = extract_normal_id(str(nf))
    if nid:
        file_mapping.append({
            'sample_id': nid,
            'group': 'normal',
            'file_path': str(nf)
        })

file_df = pd.DataFrame(file_mapping)
print(f"\n매칭된 샘플 수: {len(file_df)}")
print(f"  - 환자: {len(file_df[file_df['group']=='patient'])}")
print(f"  - 대조군: {len(file_df[file_df['group']=='normal'])}")

# 샘플 파일 저장
file_df.to_csv(OUTPUT_DIR / "sample_file_mapping.csv", index=False)
print(f"\n파일 매핑 정보 저장: {OUTPUT_DIR / 'sample_file_mapping.csv'}")

print("\n환자 ID 샘플:")
print(file_df[file_df['group']=='patient']['sample_id'].head(10).tolist())
print("\n대조군 ID 샘플:")
print(file_df[file_df['group']=='normal']['sample_id'].head(10).tolist())

# 3. 임상 데이터 정리
print("\n[3] 임상 데이터 정리 중...")

# 컬럼명 확인 및 표준화
print("\n원본 컬럼명:")
for i, col in enumerate(clinical_df.columns):
    print(f"  {i}: '{col}'")

# 임상 데이터 저장
clinical_df.to_csv(OUTPUT_DIR / "clinical_data_raw.csv", index=False)
print(f"\n원본 임상 데이터 저장: {OUTPUT_DIR / 'clinical_data_raw.csv'}")

print("\n" + "=" * 80)
print("1단계 완료: 데이터 로드 및 매칭")
print("=" * 80)
print(f"\n다음 단계:")
print("  1. 임상 데이터의 컬럼명을 확인하고 환자번호, 나이, 성별, 병기, CA19-9 컬럼 식별")
print("  2. DSS 분석 실행하여 유의미한 메틸화 영역 찾기")
print("  3. 메틸화 특징과 임상 데이터 통합")
print("  4. 예측 모델 구축")
