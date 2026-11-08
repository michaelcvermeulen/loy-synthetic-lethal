#!/usr/bin/env python3
"""
Analyze unmapped reads from HT1080 D0 samples to understand low mapping rate
"""

import gzip
import re
from collections import Counter, defaultdict
import sys
import argparse
from pathlib import Path

def load_brunello_library(library_file):
    """Load the Brunello library sgRNA sequences"""
    sgrnas = set()
    with open(library_file) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sgrna_seq = parts[1].upper()
                sgrnas.add(sgrna_seq)
    return sgrnas

def reverse_complement(seq):
    """Generate reverse complement of a sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in seq[::-1])

def extract_sgrna_at_position(seq, start, length=20):
    """Extract potential sgRNA at specified position"""
    if start + length <= len(seq):
        return seq[start:start+length].upper()
    return None

def analyze_reads(fastq_file, library_sgrnas, num_reads=100000, trim_5=22):
    """Analyze reads from FASTQ file"""

    # Expected sequences
    u6_promoter = "GTGGAAAGGACGAAACACCG"
    scaffold = "GTTTTAGAGCTAGAAATA"

    stats = {
        'total_reads': 0,
        'has_u6': 0,
        'has_scaffold': 0,
        'has_both': 0,
        'mapped_at_expected': 0,
        'mapped_at_other_positions': 0,
        'unmapped_total': 0,
        'perfect_structure_unmapped': 0,
        'one_mismatch': 0,
        'two_mismatches': 0
    }

    unmapped_sgrnas = Counter()
    mapped_sgrnas = Counter()
    position_mapped = defaultdict(int)

    with gzip.open(fastq_file, 'rt') as f:
        line_count = 0
        for line in f:
            if line_count % 4 == 1:  # Sequence line
                stats['total_reads'] += 1

                if stats['total_reads'] > num_reads:
                    break

                seq = line.strip()

                # Check for expected structural elements
                if u6_promoter in seq:
                    stats['has_u6'] += 1
                if scaffold in seq:
                    stats['has_scaffold'] += 1
                if u6_promoter in seq and scaffold in seq:
                    stats['has_both'] += 1

                # Extract sgRNA at expected position (after trim)
                if trim_5 > 0 and trim_5 < len(seq):
                    trimmed_seq = seq[trim_5:]
                    expected_sgrna = extract_sgrna_at_position(trimmed_seq, 0, 20)
                else:
                    expected_sgrna = extract_sgrna_at_position(seq, trim_5, 20)

                # Check if maps at expected position
                if expected_sgrna and expected_sgrna in library_sgrnas:
                    stats['mapped_at_expected'] += 1
                    mapped_sgrnas[expected_sgrna] += 1
                else:
                    # Try to find sgRNA at other positions
                    found = False
                    for pos in range(0, min(len(seq)-19, 50)):
                        potential_sgrna = extract_sgrna_at_position(seq, pos, 20)
                        if potential_sgrna and potential_sgrna in library_sgrnas:
                            stats['mapped_at_other_positions'] += 1
                            position_mapped[pos] += 1
                            found = True
                            break

                    if not found:
                        stats['unmapped_total'] += 1
                        if expected_sgrna:
                            unmapped_sgrnas[expected_sgrna] += 1

                            # Check if has perfect structure but unmapped sgRNA
                            if u6_promoter in seq and scaffold in seq:
                                stats['perfect_structure_unmapped'] += 1

                            # Check for mismatches
                            for lib_sgrna in library_sgrnas:
                                mismatches = sum(1 for a, b in zip(expected_sgrna, lib_sgrna) if a != b)
                                if mismatches == 1:
                                    stats['one_mismatch'] += 1
                                    break
                                elif mismatches == 2:
                                    stats['two_mismatches'] += 1
                                    break

                if stats['total_reads'] % 10000 == 0:
                    print(f"Processed {stats['total_reads']} reads...", file=sys.stderr)

            line_count += 1

    return stats, unmapped_sgrnas, mapped_sgrnas, position_mapped

def main():
    parser = argparse.ArgumentParser(description='Analyze unmapped reads from D0 samples')
    parser.add_argument('fastq', help='FASTQ file to analyze')
    parser.add_argument('library', help='Brunello library file')
    parser.add_argument('--num-reads', type=int, default=100000, help='Number of reads to analyze')
    parser.add_argument('--trim-5', type=int, default=22, help='5\' trimming length')
    parser.add_argument('--output', help='Output file for detailed results')

    args = parser.parse_args()

    print("Loading Brunello library...")
    library_sgrnas = load_brunello_library(args.library)
    print(f"Loaded {len(library_sgrnas)} sgRNAs from library")

    print(f"\nAnalyzing {args.num_reads} reads from {args.fastq}...")
    print(f"Using trim-5={args.trim_5}")

    stats, unmapped_sgrnas, mapped_sgrnas, position_mapped = analyze_reads(
        args.fastq, library_sgrnas, args.num_reads, args.trim_5
    )

    # Print summary statistics
    print("\n=== MAPPING STATISTICS ===")
    print(f"Total reads analyzed: {stats['total_reads']:,}")
    print(f"Reads with U6 promoter: {stats['has_u6']:,} ({100*stats['has_u6']/stats['total_reads']:.1f}%)")
    print(f"Reads with scaffold: {stats['has_scaffold']:,} ({100*stats['has_scaffold']/stats['total_reads']:.1f}%)")
    print(f"Reads with both U6 and scaffold: {stats['has_both']:,} ({100*stats['has_both']/stats['total_reads']:.1f}%)")

    print(f"\nMapped at expected position (trim-5={args.trim_5}): {stats['mapped_at_expected']:,} ({100*stats['mapped_at_expected']/stats['total_reads']:.1f}%)")
    print(f"Mapped at other positions: {stats['mapped_at_other_positions']:,} ({100*stats['mapped_at_other_positions']/stats['total_reads']:.1f}%)")
    print(f"Total unmapped: {stats['unmapped_total']:,} ({100*stats['unmapped_total']/stats['total_reads']:.1f}%)")

    print("\n=== UNMAPPED READ ANALYSIS ===")
    print(f"Unmapped with perfect structure (U6+scaffold): {stats['perfect_structure_unmapped']:,}")
    print(f"Unmapped that are 1 mismatch away: {stats['one_mismatch']:,}")
    print(f"Unmapped that are 2 mismatches away: {stats['two_mismatches']:,}")

    if position_mapped:
        print("\n=== ALTERNATIVE MAPPING POSITIONS ===")
        for pos, count in sorted(position_mapped.items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"Position {pos}: {count:,} reads")

    print("\n=== TOP UNMAPPED SEQUENCES ===")
    for sgrna, count in unmapped_sgrnas.most_common(20):
        print(f"{sgrna}: {count:,} reads")

    if args.output:
        with open(args.output, 'w') as f:
            f.write("# D0 Unmapped Read Analysis\n\n")
            f.write(f"File: {args.fastq}\n")
            f.write(f"Trim-5: {args.trim_5}\n")
            f.write(f"Reads analyzed: {stats['total_reads']}\n\n")

            f.write("## Statistics\n")
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")

            f.write("\n## Top 100 Unmapped sgRNAs\n")
            for sgrna, count in unmapped_sgrnas.most_common(100):
                f.write(f"{sgrna}\t{count}\n")

if __name__ == "__main__":
    main()