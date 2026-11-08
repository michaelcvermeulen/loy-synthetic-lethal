#!/usr/bin/env python3
"""
Align unmapped D0 sgRNA sequences to various CRISPR libraries
to identify contamination source
"""

import csv
import os
from collections import Counter, defaultdict
from pathlib import Path

def load_unmapped_sequences(csv_file):
    """Load unmapped sgRNA sequences from CSV file"""
    sequences = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row['extracted_20bp_sgRNA'].strip().upper()
            if len(seq) == 20 and 'N' not in seq:
                sequences.append(seq)
    return sequences

def load_brunello_library(library_file):
    """Load Brunello library sgRNA sequences"""
    sequences = {}
    with open(library_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                sgrna_id = parts[0]
                sgrna_seq = parts[1].upper()
                gene = parts[2]
                sequences[sgrna_seq] = (sgrna_id, gene)
    return sequences

def load_avana_library(library_file):
    """Load Avana library from CRISPRcleanR format"""
    sequences = {}
    with open(library_file, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num == 0 and 'GENE' in line.upper():
                continue  # Skip header if present

            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gene = parts[0]
                sgrna_id = parts[1] if len(parts) > 1 else f"Avana_{line_num}"
                sgrna_seq = parts[2].upper() if len(parts) > 2 else parts[1].upper()

                # Handle different possible formats
                if len(sgrna_seq) == 20 and 'N' not in sgrna_seq:
                    sequences[sgrna_seq] = (sgrna_id, gene)
                elif len(parts) > 3:
                    # Try other columns
                    for part in parts[1:]:
                        if len(part) == 20 and all(c in 'ACGT' for c in part.upper()):
                            sequences[part.upper()] = (sgrna_id, gene)
                            break
    return sequences

def load_tkov3_library(library_file):
    """Load TKOv3 library from CRISPRcleanR format"""
    sequences = {}
    with open(library_file, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num == 0 and any(header in line.upper() for header in ['GENE', 'SGRNA', 'REAGENT']):
                continue  # Skip header

            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Format could be: gene\tsgRNA or sgRNA_id\tsequence
                if len(parts[0]) == 20 and all(c in 'ACGT' for c in parts[0].upper()):
                    # First column is sequence
                    sgrna_seq = parts[0].upper()
                    gene = parts[1] if len(parts) > 1 else "Unknown"
                    sgrna_id = f"TKO_{line_num}"
                elif len(parts) > 1 and len(parts[1]) == 20 and all(c in 'ACGT' for c in parts[1].upper()):
                    # Second column is sequence
                    sgrna_seq = parts[1].upper()
                    gene = parts[0]
                    sgrna_id = f"TKO_{line_num}"
                else:
                    # Try to find 20bp sequence in any column
                    for i, part in enumerate(parts):
                        if len(part) == 20 and all(c in 'ACGT' for c in part.upper()):
                            sgrna_seq = part.upper()
                            gene = parts[0] if i > 0 else "Unknown"
                            sgrna_id = f"TKO_{line_num}"
                            break
                    else:
                        continue

                sequences[sgrna_seq] = (sgrna_id, gene)
    return sequences

def compare_sequences(unmapped_seqs, library_seqs, library_name):
    """Compare unmapped sequences to a library"""
    matches = []
    unmapped_counter = Counter(unmapped_seqs)

    for seq, count in unmapped_counter.items():
        if seq in library_seqs:
            sgrna_id, gene = library_seqs[seq]
            matches.append({
                'sequence': seq,
                'count': count,
                'library': library_name,
                'sgrna_id': sgrna_id,
                'gene': gene
            })

    return matches

def reverse_complement(seq):
    """Generate reverse complement of a sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in seq[::-1])

def main():
    print("=" * 70)
    print("CRISPR LIBRARY CONTAMINATION ANALYSIS")
    print("=" * 70)

    # Define file paths
    base_dir = Path("/global/home/hpc5173/projects/CRISPR_screen")
    unmapped_ctrl = base_dir / "analysis/ht1080_screen_2/D0_CTRL_unmapped_sequences.csv"
    unmapped_loy = base_dir / "analysis/ht1080_screen_2/D0_LOY_unmapped_sequences.csv"

    brunello_lib = base_dir / "data/brunello/broadgpp-brunello-library-corrected.txt"
    avana_lib = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/Avana_library_forCRISPRcleanR_hg38.txt"
    tkov3_lib = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/TKOv3_library_forCRISPRcleanR_sgRNAtag.txt"

    # Load unmapped sequences
    print("\nLoading unmapped sequences...")
    ctrl_seqs = load_unmapped_sequences(unmapped_ctrl)
    loy_seqs = load_unmapped_sequences(unmapped_loy)
    all_unmapped = ctrl_seqs + loy_seqs

    print(f"  D0-CTRL: {len(ctrl_seqs)} sequences")
    print(f"  D0-LOY: {len(loy_seqs)} sequences")
    print(f"  Total unique: {len(set(all_unmapped))} sequences")

    # Load libraries
    print("\nLoading CRISPR libraries...")

    print("  Loading Brunello...")
    brunello_seqs = load_brunello_library(brunello_lib)
    print(f"    Loaded {len(brunello_seqs)} Brunello sgRNAs")

    print("  Loading Avana...")
    if avana_lib.exists():
        avana_seqs = load_avana_library(avana_lib)
        print(f"    Loaded {len(avana_seqs)} Avana sgRNAs")
    else:
        print(f"    Avana library not found at {avana_lib}")
        avana_seqs = {}

    print("  Loading TKOv3...")
    if tkov3_lib.exists():
        tkov3_seqs = load_tkov3_library(tkov3_lib)
        print(f"    Loaded {len(tkov3_seqs)} TKOv3 sgRNAs")
    else:
        print(f"    TKOv3 library not found at {tkov3_lib}")
        tkov3_seqs = {}

    # Compare sequences
    print("\n" + "=" * 70)
    print("MATCHING UNMAPPED SEQUENCES TO LIBRARIES")
    print("=" * 70)

    # Check Brunello (sanity check - should be minimal)
    brunello_matches = compare_sequences(all_unmapped, brunello_seqs, "Brunello")
    print(f"\nBrunello matches: {len(brunello_matches)} sequences")
    print("  (These shouldn't exist - indicates mapping issue)")

    # Check Avana
    avana_matches = compare_sequences(all_unmapped, avana_seqs, "Avana")
    print(f"\nAvana matches: {len(avana_matches)} sequences")

    # Check TKOv3
    tkov3_matches = compare_sequences(all_unmapped, tkov3_seqs, "TKOv3")
    print(f"TKOv3 matches: {len(tkov3_matches)} sequences")

    # Check reverse complements
    print("\nChecking reverse complements...")
    all_unmapped_rc = [reverse_complement(seq) for seq in all_unmapped]
    brunello_rc_matches = compare_sequences(all_unmapped_rc, brunello_seqs, "Brunello_RC")
    avana_rc_matches = compare_sequences(all_unmapped_rc, avana_seqs, "Avana_RC")
    tkov3_rc_matches = compare_sequences(all_unmapped_rc, tkov3_seqs, "TKOv3_RC")

    print(f"  Brunello RC matches: {len(brunello_rc_matches)}")
    print(f"  Avana RC matches: {len(avana_rc_matches)}")
    print(f"  TKOv3 RC matches: {len(tkov3_rc_matches)}")

    # Calculate statistics
    total_unmapped = len(set(all_unmapped))
    matched_seqs = set()

    for matches in [brunello_matches, avana_matches, tkov3_matches]:
        matched_seqs.update([m['sequence'] for m in matches])

    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    print(f"\nTotal unique unmapped sequences: {total_unmapped}")
    print(f"Sequences matched to any library: {len(matched_seqs)} ({100*len(matched_seqs)/total_unmapped:.1f}%)")
    print(f"Sequences with no match: {total_unmapped - len(matched_seqs)} ({100*(total_unmapped-len(matched_seqs))/total_unmapped:.1f}%)")

    # Output results
    output_dir = base_dir / "analysis/ht1080_screen_2/library_comparison"
    output_dir.mkdir(exist_ok=True)

    # Save detailed matches
    if avana_matches:
        with open(output_dir / "avana_matches.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['sequence', 'count', 'gene', 'sgrna_id'])
            writer.writeheader()
            writer.writerows(sorted(avana_matches, key=lambda x: x['count'], reverse=True))
        print(f"\nSaved Avana matches to avana_matches.csv")

    if tkov3_matches:
        with open(output_dir / "tkov3_matches.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['sequence', 'count', 'gene', 'sgrna_id'])
            writer.writeheader()
            writer.writerows(sorted(tkov3_matches, key=lambda x: x['count'], reverse=True))
        print(f"Saved TKOv3 matches to tkov3_matches.csv")

    # Print top matches
    if avana_matches:
        print("\n" + "=" * 70)
        print("TOP 20 AVANA MATCHES")
        print("=" * 70)
        sorted_avana = sorted(avana_matches, key=lambda x: x['count'], reverse=True)
        for i, match in enumerate(sorted_avana[:20], 1):
            print(f"{i:2d}. {match['sequence']}: {match['count']:3d} reads -> {match['gene']}")

    if tkov3_matches:
        print("\n" + "=" * 70)
        print("TOP 20 TKOv3 MATCHES")
        print("=" * 70)
        sorted_tkov3 = sorted(tkov3_matches, key=lambda x: x['count'], reverse=True)
        for i, match in enumerate(sorted_tkov3[:20], 1):
            print(f"{i:2d}. {match['sequence']}: {match['count']:3d} reads -> {match['gene']}")

    # Generate summary report
    with open(output_dir / "contamination_summary.txt", 'w') as f:
        f.write("CONTAMINATION SOURCE ANALYSIS\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total unique unmapped sequences: {total_unmapped}\n")
        f.write(f"Brunello matches: {len(brunello_matches)} ({100*len(brunello_matches)/total_unmapped:.1f}%)\n")
        f.write(f"Avana matches: {len(avana_matches)} ({100*len(avana_matches)/total_unmapped:.1f}%)\n")
        f.write(f"TKOv3 matches: {len(tkov3_matches)} ({100*len(tkov3_matches)/total_unmapped:.1f}%)\n")
        f.write(f"No matches: {total_unmapped - len(matched_seqs)} ({100*(total_unmapped-len(matched_seqs))/total_unmapped:.1f}%)\n")

        if len(avana_matches) > len(tkov3_matches) * 2:
            f.write("\nCONCLUSION: Contamination appears to be from AVANA library\n")
        elif len(tkov3_matches) > len(avana_matches) * 2:
            f.write("\nCONCLUSION: Contamination appears to be from TKOv3 library\n")
        elif len(matched_seqs) < total_unmapped * 0.1:
            f.write("\nCONCLUSION: Contamination source not identified in tested libraries\n")
            f.write("Consider checking GeCKO v2 or other CRISPR libraries\n")
        else:
            f.write("\nCONCLUSION: Mixed contamination from multiple sources\n")

    print(f"\nSummary saved to contamination_summary.txt")

if __name__ == "__main__":
    main()