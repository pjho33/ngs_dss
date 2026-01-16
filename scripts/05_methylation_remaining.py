# scripts/05_methylation_remaining.py
# 남은 샘플들만 methylation extraction 진행
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def run_methylation_extractor_remaining():
    paths, params = load_paths_and_params()

    output_root = Path(paths.get("output_dir", "results"))
    bam_dir = output_root / "bam"
    meth_dir = output_root / "methylation"

    threads = params.get("methylation", {}).get("threads", 8)

    # 이미 완료된 샘플 목록
    done_samples = set()
    for cov in meth_dir.glob("*.bismark.cov.gz"):
        # 파일명에서 샘플명 추출
        sample_name = cov.name.replace(".bismark.cov.gz", "")
        done_samples.add(sample_name)

    # 모든 deduplicated BAM 파일
    all_bams = sorted(bam_dir.glob("*.deduplicated.bam"))
    
    # 아직 처리 안된 BAM 파일만 필터링
    remaining_bams = []
    for bam in all_bams:
        sample_name = bam.stem  # .bam 제거
        if sample_name not in done_samples:
            remaining_bams.append(bam)

    print(f"[METH] Total BAMs: {len(all_bams)}, Done: {len(done_samples)}, Remaining: {len(remaining_bams)}")

    if not remaining_bams:
        print("[METH] All samples already processed!")
        return

    for i, bam in enumerate(remaining_bams, 1):
        print(f"[METH] Processing {i}/{len(remaining_bams)}: {bam.name}")
        cmd = [
            "conda", "run", "-n", "bismark_env",
            "bismark_methylation_extractor",
            "--bedGraph",
            "--gzip",
            "--multicore", str(threads),
            "--comprehensive",
            "--cytosine_report",
            "--genome_folder", str(Path(paths["reference_genome"]).parent),
            "-o", str(meth_dir),
            str(bam),
        ]
        run_cmd(cmd, log_name=f"05_methylation_extractor_{bam.stem}.log")

    print(f"[METH] Methylation extraction finished. Results in: {meth_dir}")


if __name__ == "__main__":
    run_methylation_extractor_remaining()
