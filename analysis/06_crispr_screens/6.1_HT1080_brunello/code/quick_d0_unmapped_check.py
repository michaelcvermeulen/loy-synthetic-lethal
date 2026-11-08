#!/usr/bin/env python3
"""
Quick analysis of unmapped reads from D0 samples
"""

import gzip
from collections import Counter

# Expected sequences
U6_PROMOTER = "GTGGAAAGGACGAAACACCG"
SCAFFOLD = "GTTTTAGAGCTAGAAATA"

def analyze_d0_reads(fastq_file, num_reads=10000):
    """Quick analysis of D0 reads"""

    stats = {
        'total': 0,
        'has_u6': 0,
        'has_scaffold': 0,
        'has_both': 0,
        'u6_at_pos1': 0,
        'u6_at_pos2': 0,
        'scaffold_after_20bp': 0
    }

    sgrna_sequences = Counter()

    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                seq = line.strip()
                stats['total'] += 1

                if stats['total'] > num_reads:
                    break

                # Check for structural elements
                if U6_PROMOTER in seq:
                    stats['has_u6'] += 1
                    u6_pos = seq.find(U6_PROMOTER)
                    if u6_pos == 1:
                        stats['u6_at_pos1'] += 1
                    elif u6_pos == 2:
                        stats['u6_at_pos2'] += 1

                if SCAFFOLD in seq:
                    stats['has_scaffold'] += 1
                    scaffold_pos = seq.find(SCAFFOLD)
                    # Check if scaffold is ~20bp after U6
                    if U6_PROMOTER in seq:
                        u6_pos = seq.find(U6_PROMOTER)
                        if 38 <= (scaffold_pos - u6_pos) <= 42:
                            stats['scaffold_after_20bp'] += 1

                if U6_PROMOTER in seq and SCAFFOLD in seq:
                    stats['has_both'] += 1
                    # Extract potential sgRNA
                    u6_pos = seq.find(U6_PROMOTER)
                    scaffold_pos = seq.find(SCAFFOLD)
                    if u6_pos >= 0 and scaffold_pos > u6_pos + 20:
                        sgrna_start = u6_pos + len(U6_PROMOTER)
                        sgrna_seq = seq[sgrna_start:sgrna_start+20]
                        if len(sgrna_seq) == 20:
                            sgrna_sequences[sgrna_seq] += 1

    return stats, sgrna_sequences

# Analyze D0-CTRL
print("Analyzing D0-CTRL reads...")
d0_ctrl_file = "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/D0-CTRL_R1_001.fastq.gz"
stats, sgrnas = analyze_d0_reads(d0_ctrl_file, 10000)

print(f"\n=== D0-CTRL ANALYSIS (first 10,000 reads) ===")
print(f"Total reads: {stats['total']:,}")
print(f"Has U6 promoter: {stats['has_u6']:,} ({100*stats['has_u6']/stats['total']:.1f}%)")
print(f"Has scaffold: {stats['has_scaffold']:,} ({100*stats['has_scaffold']/stats['total']:.1f}%)")
print(f"Has both U6 and scaffold: {stats['has_both']:,} ({100*stats['has_both']/stats['total']:.1f}%)")
print(f"\n=== STRUCTURAL POSITIONS ===")
print(f"U6 at position 1: {stats['u6_at_pos1']:,} ({100*stats['u6_at_pos1']/stats['total']:.1f}%)")
print(f"U6 at position 2: {stats['u6_at_pos2']:,} ({100*stats['u6_at_pos2']/stats['total']:.1f}%)")
print(f"Scaffold ~20bp after U6: {stats['scaffold_after_20bp']:,} ({100*stats['scaffold_after_20bp']/stats['total']:.1f}%)")

print(f"\n=== TOP 10 EXTRACTED sgRNAs ===")
for sgrna, count in sgrnas.most_common(10):
    print(f"{sgrna}: {count} reads")

# Compare with POST sample
print("\n" + "="*60)
print("Analyzing POST-CTRL-REP1 reads for comparison...")
post_file = "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/CTRL-REP1_R1_001.fastq.gz"
post_stats, post_sgrnas = analyze_d0_reads(post_file, 10000)

print(f"\n=== POST-CTRL-REP1 ANALYSIS (first 10,000 reads) ===")
print(f"Total reads: {post_stats['total']:,}")
print(f"Has U6 promoter: {post_stats['has_u6']:,} ({100*post_stats['has_u6']/post_stats['total']:.1f}%)")
print(f"Has scaffold: {post_stats['has_scaffold']:,} ({100*post_stats['has_scaffold']/post_stats['total']:.1f}%)")
print(f"Has both U6 and scaffold: {post_stats['has_both']:,} ({100*post_stats['has_both']/post_stats['total']:.1f}%)")
print(f"\n=== STRUCTURAL POSITIONS ===")
print(f"U6 at position 1: {post_stats['u6_at_pos1']:,} ({100*post_stats['u6_at_pos1']/post_stats['total']:.1f}%)")
print(f"U6 at position 2: {post_stats['u6_at_pos2']:,} ({100*post_stats['u6_at_pos2']/post_stats['total']:.1f}%)")
print(f"Scaffold ~20bp after U6: {post_stats['scaffold_after_20bp']:,} ({100*post_stats['scaffold_after_20bp']/post_stats['total']:.1f}%)")