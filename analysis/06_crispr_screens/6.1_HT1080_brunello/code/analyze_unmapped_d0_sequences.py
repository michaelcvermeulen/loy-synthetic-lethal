#!/usr/bin/env python3
"""
Extract and analyze sgRNA sequences from D0 reads to understand low mapping
"""

import gzip
from collections import Counter
import random

def extract_all_sgrnas(fastq_file, num_reads=100000):
    """Extract sgRNA sequences from reads regardless of mapping"""

    U6_PROMOTER = "GTGGAAAGGACGAAACACCG"
    SCAFFOLD = "GTTTTAGAGCTAGAAATA"

    all_sgrnas = Counter()
    no_structure = 0
    total = 0

    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                seq = line.strip()
                total += 1

                if total > num_reads:
                    break

                # Find U6 and scaffold positions
                u6_pos = seq.find(U6_PROMOTER)
                scaffold_pos = seq.find(SCAFFOLD)

                if u6_pos >= 0 and scaffold_pos > u6_pos:
                    # Extract 20bp between U6 and scaffold
                    sgrna_start = u6_pos + len(U6_PROMOTER)
                    sgrna_end = scaffold_pos

                    if sgrna_end - sgrna_start >= 20:
                        sgrna = seq[sgrna_start:sgrna_start + 20]
                        if len(sgrna) == 20 and 'N' not in sgrna:
                            all_sgrnas[sgrna] += 1
                else:
                    no_structure += 1

    return all_sgrnas, total, no_structure

# Load Brunello library
brunello_sgrnas = set()
with open("/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt") as f:
    f.readline()  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            brunello_sgrnas.add(parts[1].upper())

print("=" * 70)
print("D0 UNMAPPED SEQUENCE ANALYSIS")
print("=" * 70)

# Analyze D0 sample
print("\nExtracting sgRNAs from D0-CTRL...")
d0_sgrnas, d0_total, d0_no_struct = extract_all_sgrnas(
    "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/D0-CTRL_R1_001.fastq.gz",
    100000
)

# Analyze POST sample for comparison
print("Extracting sgRNAs from POST-CTRL-REP1...")
post_sgrnas, post_total, post_no_struct = extract_all_sgrnas(
    "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/CTRL-REP1_R1_001.fastq.gz",
    100000
)

print("\n" + "=" * 70)
print("EXTRACTION STATISTICS")
print("=" * 70)

print(f"\nD0-CTRL:")
print(f"  Total reads: {d0_total:,}")
print(f"  Reads with proper structure: {d0_total - d0_no_struct:,} ({100*(d0_total-d0_no_struct)/d0_total:.1f}%)")
print(f"  Unique sgRNA sequences found: {len(d0_sgrnas):,}")

print(f"\nPOST-CTRL:")
print(f"  Total reads: {post_total:,}")
print(f"  Reads with proper structure: {post_total - post_no_struct:,} ({100*(post_total-post_no_struct)/post_total:.1f}%)")
print(f"  Unique sgRNA sequences found: {len(post_sgrnas):,}")

# Check mapping to Brunello
d0_in_brunello = sum(1 for s in d0_sgrnas if s in brunello_sgrnas)
d0_not_in_brunello = sum(1 for s in d0_sgrnas if s not in brunello_sgrnas)

post_in_brunello = sum(1 for s in post_sgrnas if s in brunello_sgrnas)
post_not_in_brunello = sum(1 for s in post_sgrnas if s not in brunello_sgrnas)

print("\n" + "=" * 70)
print("BRUNELLO LIBRARY MAPPING")
print("=" * 70)

print(f"\nD0 unique sgRNAs:")
print(f"  In Brunello library: {d0_in_brunello:,} ({100*d0_in_brunello/len(d0_sgrnas):.1f}%)")
print(f"  NOT in library: {d0_not_in_brunello:,} ({100*d0_not_in_brunello/len(d0_sgrnas):.1f}%)")

print(f"\nPOST unique sgRNAs:")
print(f"  In Brunello library: {post_in_brunello:,} ({100*post_in_brunello/len(post_sgrnas):.1f}%)")
print(f"  NOT in library: {post_not_in_brunello:,} ({100*post_not_in_brunello/len(post_sgrnas):.1f}%)")

# Count total reads for mapped vs unmapped
d0_mapped_reads = sum(count for sgrna, count in d0_sgrnas.items() if sgrna in brunello_sgrnas)
d0_unmapped_reads = sum(count for sgrna, count in d0_sgrnas.items() if sgrna not in brunello_sgrnas)

post_mapped_reads = sum(count for sgrna, count in post_sgrnas.items() if sgrna in brunello_sgrnas)
post_unmapped_reads = sum(count for sgrna, count in post_sgrnas.items() if sgrna not in brunello_sgrnas)

print("\n" + "=" * 70)
print("READ-LEVEL MAPPING")
print("=" * 70)

print(f"\nD0 reads with extracted sgRNAs:")
print(f"  Mapped to Brunello: {d0_mapped_reads:,} ({100*d0_mapped_reads/(d0_mapped_reads+d0_unmapped_reads):.1f}%)")
print(f"  NOT mapped: {d0_unmapped_reads:,} ({100*d0_unmapped_reads/(d0_mapped_reads+d0_unmapped_reads):.1f}%)")

print(f"\nPOST reads with extracted sgRNAs:")
print(f"  Mapped to Brunello: {post_mapped_reads:,} ({100*post_mapped_reads/(post_mapped_reads+post_unmapped_reads):.1f}%)")
print(f"  NOT mapped: {post_unmapped_reads:,} ({100*post_unmapped_reads/(post_mapped_reads+post_unmapped_reads):.1f}%)")

# Show top unmapped sequences
print("\n" + "=" * 70)
print("TOP UNMAPPED sgRNA SEQUENCES IN D0")
print("=" * 70)

unmapped_d0 = [(s, c) for s, c in d0_sgrnas.items() if s not in brunello_sgrnas]
unmapped_d0.sort(key=lambda x: x[1], reverse=True)

print("\nTop 20 non-Brunello sgRNAs in D0:")
for i, (sgrna, count) in enumerate(unmapped_d0[:20], 1):
    # Check if it's one mismatch away from a Brunello guide
    close_match = None
    for brun_sgrna in random.sample(list(brunello_sgrnas), min(1000, len(brunello_sgrnas))):
        mismatches = sum(1 for a, b in zip(sgrna, brun_sgrna) if a != b)
        if mismatches == 1:
            close_match = brun_sgrna
            break

    if close_match:
        print(f"{i:2d}. {sgrna}: {count:4d} reads (1 mismatch from library)")
    else:
        print(f"{i:2d}. {sgrna}: {count:4d} reads")

# Check overlap between D0 and POST unmapped
d0_unmapped_set = set(s for s in d0_sgrnas if s not in brunello_sgrnas)
post_unmapped_set = set(s for s in post_sgrnas if s not in brunello_sgrnas)
overlap = d0_unmapped_set & post_unmapped_set

print(f"\nUnmapped sgRNAs shared between D0 and POST: {len(overlap)}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("\nThe analysis shows that D0 samples have:")
print(f"1. {100*(d0_total-d0_no_struct)/d0_total:.1f}% of reads with correct structure")
print(f"2. Only {100*d0_mapped_reads/(d0_mapped_reads+d0_unmapped_reads):.1f}% of extracted sgRNAs map to Brunello")
print(f"3. POST samples have {100*post_mapped_reads/(post_mapped_reads+post_unmapped_reads):.1f}% mapping")
print("\nThis suggests the unmapped D0 reads contain non-Brunello sgRNA sequences!")