#!/usr/bin/env python3
"""
Extract unmapped reads from D0 samples and output to CSV
"""

import gzip
import csv
import sys

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

def extract_unmapped_sequences(fastq_file, library_sgrnas, output_file, num_unmapped=1000):
    """Extract unmapped sequences and save to CSV"""

    U6_PROMOTER = "GTGGAAAGGACGAAACACCG"
    SCAFFOLD = "GTTTTAGAGCTAGAAATA"

    unmapped_found = 0
    total_processed = 0

    with gzip.open(fastq_file, 'rt') as f, open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['read_id', 'full_sequence', 'extracted_20bp_sgRNA', 'u6_position', 'scaffold_position'])

        lines = []
        for line in f:
            lines.append(line.strip())

            if len(lines) == 4:  # Complete read
                read_id = lines[0][1:]  # Remove @ symbol
                sequence = lines[1]

                total_processed += 1

                # Find U6 and scaffold positions
                u6_pos = sequence.find(U6_PROMOTER)
                scaffold_pos = sequence.find(SCAFFOLD)

                # Extract sgRNA if structure is correct
                extracted_sgrna = ""
                if u6_pos >= 0 and scaffold_pos > u6_pos:
                    sgrna_start = u6_pos + len(U6_PROMOTER)
                    if scaffold_pos - sgrna_start >= 20:
                        extracted_sgrna = sequence[sgrna_start:sgrna_start + 20]

                # Alternative extraction with trim-5=22 (MAGeCK default)
                if not extracted_sgrna and len(sequence) >= 42:
                    trimmed_seq = sequence[22:]
                    alt_sgrna = trimmed_seq[:20] if len(trimmed_seq) >= 20 else ""

                    # Check if this alternative extraction is in library
                    if alt_sgrna and alt_sgrna not in library_sgrnas and 'N' not in alt_sgrna:
                        extracted_sgrna = alt_sgrna

                # Check if unmapped
                if extracted_sgrna and extracted_sgrna not in library_sgrnas and 'N' not in extracted_sgrna:
                    writer.writerow([read_id, sequence, extracted_sgrna, u6_pos, scaffold_pos])
                    unmapped_found += 1

                    if unmapped_found >= num_unmapped:
                        break

                # Clear for next read
                lines = []

                # Progress update
                if total_processed % 10000 == 0:
                    print(f"  Processed {total_processed:,} reads, found {unmapped_found:,} unmapped...", end='\r')

        print(f"\n  Final: Processed {total_processed:,} reads, found {unmapped_found:,} unmapped sequences")

    return unmapped_found, total_processed

def main():
    # Load Brunello library
    library_file = "/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt"
    print("Loading Brunello library...")
    library_sgrnas = load_brunello_library(library_file)
    print(f"Loaded {len(library_sgrnas)} sgRNAs from library\n")

    # Process D0-CTRL
    print("Processing D0-CTRL...")
    d0_ctrl_file = "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/D0-CTRL_R1_001.fastq.gz"
    d0_ctrl_output = "/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/D0_CTRL_unmapped_sequences.csv"

    ctrl_unmapped, ctrl_total = extract_unmapped_sequences(
        d0_ctrl_file,
        library_sgrnas,
        d0_ctrl_output,
        num_unmapped=1000
    )

    print(f"D0-CTRL: Saved {ctrl_unmapped} unmapped sequences to {d0_ctrl_output}")
    print(f"         Unmapped rate: {100*ctrl_unmapped/ctrl_total:.1f}% in first {ctrl_total:,} reads\n")

    # Process D0-LOY
    print("Processing D0-LOY...")
    d0_loy_file = "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/D0-LOY_R1_001.fastq.gz"
    d0_loy_output = "/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/D0_LOY_unmapped_sequences.csv"

    loy_unmapped, loy_total = extract_unmapped_sequences(
        d0_loy_file,
        library_sgrnas,
        d0_loy_output,
        num_unmapped=1000
    )

    print(f"D0-LOY: Saved {loy_unmapped} unmapped sequences to {d0_loy_output}")
    print(f"        Unmapped rate: {100*loy_unmapped/loy_total:.1f}% in first {loy_total:,} reads\n")

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nOutput files created:")
    print(f"1. {d0_ctrl_output}")
    print(f"   - Contains 1000 unmapped reads from D0-CTRL")
    print(f"   - Columns: read_id, full_sequence, extracted_20bp_sgRNA, u6_position, scaffold_position")
    print(f"\n2. {d0_loy_output}")
    print(f"   - Contains 1000 unmapped reads from D0-LOY")
    print(f"   - Same column structure")
    print(f"\nYou can now check these sequences against other CRISPR libraries or databases.")

if __name__ == "__main__":
    main()