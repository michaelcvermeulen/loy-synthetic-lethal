#!/usr/bin/env python3
"""
Test different trim positions to optimize D0 mapping
"""

import gzip
import subprocess
import os
from pathlib import Path

def create_test_fastq(input_file, output_file, num_reads=1000000):
    """Create a smaller test FASTQ for faster testing"""
    with gzip.open(input_file, 'rt') as inf, gzip.open(output_file, 'wt') as outf:
        line_count = 0
        for line in inf:
            outf.write(line)
            line_count += 1
            if line_count >= num_reads * 4:  # 4 lines per read
                break
    print(f"Created test file with {num_reads} reads: {output_file}")

def test_trim_parameter(fastq_file, library_file, trim_5, output_prefix):
    """Test MAGeCK count with specific trim parameter"""

    cmd = [
        "mageck", "count",
        "-l", library_file,
        "--sample-label", f"trim_{trim_5}",
        "--fastq", fastq_file,
        "--trim-5", str(trim_5),
        "-n", output_prefix,
        "--norm-method", "none"  # Skip normalization for speed
    ]

    print(f"\nTesting trim-5={trim_5}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse the log for mapping statistics
    for line in result.stderr.split('\n'):
        if 'Total:' in line or 'Mapped:' in line:
            print(f"  {line.strip()}")

    # Get mapping stats from log file
    log_file = f"{output_prefix}.log"
    if os.path.exists(log_file):
        with open(log_file) as f:
            for line in f:
                if 'mappedreads' in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        mapped = int(parts[1])
                        return mapped
    return 0

# Create test directory
test_dir = Path("/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/trim_optimization")
test_dir.mkdir(exist_ok=True)

# Create test FASTQ files
print("Creating test FASTQ files...")
d0_test = test_dir / "d0_ctrl_test.fastq.gz"
post_test = test_dir / "post_ctrl_test.fastq.gz"

create_test_fastq(
    "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/D0-CTRL_R1_001.fastq.gz",
    d0_test,
    100000
)

create_test_fastq(
    "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/CTRL-REP1_R1_001.fastq.gz",
    post_test,
    100000
)

# Test different trim parameters
library = "/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt"

trim_values = [0, 1, 2, 10, 15, 20, 21, 22, 23, 24, 25, 30]

print("\n" + "="*60)
print("TESTING D0-CTRL WITH DIFFERENT TRIM VALUES")
print("="*60)

d0_results = {}
for trim in trim_values:
    output_prefix = test_dir / f"d0_trim_{trim}"
    mapped = test_trim_parameter(str(d0_test), library, trim, str(output_prefix))
    d0_results[trim] = mapped
    print(f"  Mapped reads: {mapped:,} ({100*mapped/100000:.1f}%)")

print("\n" + "="*60)
print("TESTING POST-CTRL WITH DIFFERENT TRIM VALUES")
print("="*60)

post_results = {}
for trim in trim_values:
    output_prefix = test_dir / f"post_trim_{trim}"
    mapped = test_trim_parameter(str(post_test), library, trim, str(output_prefix))
    post_results[trim] = mapped
    print(f"  Mapped reads: {mapped:,} ({100*mapped/100000:.1f}%)")

print("\n" + "="*60)
print("SUMMARY COMPARISON")
print("="*60)
print(f"{'Trim-5':<10} {'D0 Mapped %':<15} {'POST Mapped %':<15} {'Ratio'}")
print("-" * 50)
for trim in trim_values:
    d0_pct = 100 * d0_results.get(trim, 0) / 100000
    post_pct = 100 * post_results.get(trim, 0) / 100000
    ratio = post_pct / d0_pct if d0_pct > 0 else 0
    print(f"{trim:<10} {d0_pct:<15.1f} {post_pct:<15.1f} {ratio:.2f}x")