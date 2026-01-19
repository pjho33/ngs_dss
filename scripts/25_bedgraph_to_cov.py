#!/usr/bin/env python3
"""
bedGraph 파일을 Bismark .cov 형식으로 변환
bedGraph: chr start end methylation_percentage
cov: chr start end methylation_percentage count_methylated count_unmethylated
"""

import gzip
import os
from pathlib import Path

# Input and output directories
input_dir = Path("/home/pjho3/바탕화면/normalBedgraph")
output_dir = Path("/home/pjho3/바탕화면/normal_cov")
output_dir.mkdir(exist_ok=True)

print("=" * 80)
print("bedGraph to .cov 변환")
print("=" * 80)

# Get all bedGraph files (exclude hidden files starting with ._)
bedgraph_files = [f for f in sorted(input_dir.glob("*.bedGraph.gz")) 
                  if not f.name.startswith('._')]
print(f"\n발견된 bedGraph 파일: {len(bedgraph_files)}개")

for idx, bedgraph_file in enumerate(bedgraph_files, 1):
    # Create output filename
    sample_name = bedgraph_file.name.replace(".bedGraph.gz", "")
    output_file = output_dir / f"{sample_name}.bismark.cov.gz"
    
    print(f"\n[{idx}/{len(bedgraph_files)}] 변환 중: {bedgraph_file.name}")
    
    line_count = 0
    
    with gzip.open(bedgraph_file, 'rt') as fin, gzip.open(output_file, 'wt') as fout:
        for line in fin:
            line = line.strip()
            
            # Skip header
            if line.startswith('track'):
                continue
            
            parts = line.split('\t')
            if len(parts) != 4:
                continue
            
            chrom, start, end, meth_pct = parts
            
            # bedGraph는 0-based, cov는 1-based이므로 start에 +1
            start = int(start) + 1
            end = int(end)
            meth_pct = float(meth_pct)
            
            # 메틸화 비율을 count로 변환 (가정: coverage = 10)
            # 실제 coverage 정보가 없으므로 임의로 10으로 설정
            coverage = 10
            count_methylated = int(round(meth_pct * coverage / 100))
            count_unmethylated = coverage - count_methylated
            
            # .cov 형식: chr start end methylation_percentage count_methylated count_unmethylated
            fout.write(f"{chrom}\t{start}\t{end}\t{meth_pct:.2f}\t{count_methylated}\t{count_unmethylated}\n")
            line_count += 1
    
    print(f"   완료: {line_count:,} lines -> {output_file.name}")

print("\n" + "=" * 80)
print(f"변환 완료! 출력 디렉토리: {output_dir}")
print("=" * 80)

# Summary
print(f"\n변환된 파일 목록:")
for cov_file in sorted(output_dir.glob("*.cov.gz")):
    size_mb = cov_file.stat().st_size / (1024 * 1024)
    print(f"  {cov_file.name} ({size_mb:.1f} MB)")
