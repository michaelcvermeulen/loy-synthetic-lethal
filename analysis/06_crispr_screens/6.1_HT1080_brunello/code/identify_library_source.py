#!/usr/bin/env python3
"""
Properly identify the source library of unmapped sequences by exact matching.
Includes genomic coordinates for all matched sequences.
"""

import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

def load_library_with_coordinates(file_path, lib_name):
    """Load a CRISPR library file and extract sequences with coordinates"""
    sequences = {}

    if lib_name == "Avana":
        # Two possible Avana files
        if "forCRISPRcleanR" in str(file_path):
            # This file has genomic coordinates
            df = pd.read_csv(file_path, sep='\t')
            for _, row in df.iterrows():
                seq = row['CODE'].upper().strip()
                sequences[seq] = {
                    'gene': row['GENES'],
                    'exon': row.get('EXONE', 'NA'),
                    'chromosome': str(row.get('CHRM', 'NA')),
                    'strand': row.get('STRAND', 'NA'),
                    'start': row.get('STARTpos', 'NA'),
                    'end': row.get('ENDpos', 'NA')
                }
        else:
            # MAGeCK format file
            with open(file_path, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        seq = parts[1].upper().strip()
                        gene = parts[2] if len(parts) > 2 else 'Unknown'
                        sequences[seq] = {
                            'gene': gene,
                            'chromosome': 'NA',
                            'start': 'NA',
                            'end': 'NA'
                        }

    elif lib_name == "Brunello":
        with open(file_path, 'r') as f:
            # Check if there's a header
            first_line = f.readline().strip()
            if not first_line.startswith('sgRNA'):
                # No header, process first line
                parts = first_line.split('\t')
                if len(parts) >= 3:
                    seq = parts[1].upper().strip()
                    gene = parts[2]
                    sequences[seq] = {'gene': gene, 'chromosome': 'NA', 'start': 'NA', 'end': 'NA'}

            # Process rest of file
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    seq = parts[1].upper().strip()
                    gene = parts[2]
                    sequences[seq] = {'gene': gene, 'chromosome': 'NA', 'start': 'NA', 'end': 'NA'}

    elif lib_name == "TKOv3":
        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            if 'CODE' in row:
                seq = str(row['CODE']).upper().strip()
                sequences[seq] = {
                    'gene': row.get('GENES', 'NA'),
                    'chromosome': str(row.get('CHRM', 'NA')),
                    'start': row.get('STARTpos', 'NA'),
                    'end': row.get('ENDpos', 'NA'),
                    'strand': row.get('STRAND', 'NA')
                }

    return sequences

def main():
    base_dir = Path("/global/home/hpc5173/projects/CRISPR_screen")
    output_dir = base_dir / "analysis/ht1080_screen_2/genomic_mapping"
    output_dir.mkdir(exist_ok=True)

    print("="*70)
    print("LIBRARY SOURCE IDENTIFICATION WITH GENOMIC COORDINATES")
    print("="*70)

    # Load unmapped sequences
    print("\nLoading unmapped sequences...")
    ctrl_file = base_dir / "analysis/ht1080_screen_2/D0_CTRL_unmapped_sequences.csv"
    loy_file = base_dir / "analysis/ht1080_screen_2/D0_LOY_unmapped_sequences.csv"

    ctrl_df = pd.read_csv(ctrl_file)
    loy_df = pd.read_csv(loy_file)

    # Get counts
    ctrl_sequences = ctrl_df['extracted_20bp_sgRNA'].dropna().str.upper().str.strip()
    loy_sequences = loy_df['extracted_20bp_sgRNA'].dropna().str.upper().str.strip()

    ctrl_counts = Counter(ctrl_sequences)
    loy_counts = Counter(loy_sequences)

    # Combine all unique sequences
    all_sequences = set(ctrl_counts.keys()) | set(loy_counts.keys())

    print(f"Found {len(all_sequences)} unique unmapped sequences")
    print(f"Total unmapped reads: {sum(ctrl_counts.values()) + sum(loy_counts.values())}")

    # Load all libraries
    print("\nLoading CRISPR libraries...")

    # Avana with coordinates
    avana_file1 = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/Avana_library_forCRISPRcleanR_hg38.txt"
    avana_seqs = load_library_with_coordinates(avana_file1, "Avana")
    print(f"Loaded {len(avana_seqs)} Avana sequences with coordinates")

    # Also load converted Avana
    avana_file2 = base_dir / "data/avana/Avana_library_MAGeCK.txt"
    if avana_file2.exists():
        avana_seqs2 = load_library_with_coordinates(avana_file2, "Avana")
        # Merge, keeping sequences with coordinates when available
        for seq, info in avana_seqs2.items():
            if seq not in avana_seqs:
                avana_seqs[seq] = info
        print(f"Total Avana sequences after merging: {len(avana_seqs)}")

    # Brunello
    brunello_file = base_dir / "data/brunello/broadgpp-brunello-library-corrected.txt"
    brunello_seqs = load_library_with_coordinates(brunello_file, "Brunello")
    print(f"Loaded {len(brunello_seqs)} Brunello sequences")

    # TKOv3
    tkov3_file = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/TKOv3_library_forCRISPRcleanR_sgRNAtag.txt"
    tkov3_seqs = {}
    if tkov3_file.exists():
        tkov3_seqs = load_library_with_coordinates(tkov3_file, "TKOv3")
        print(f"Loaded {len(tkov3_seqs)} TKOv3 sequences")

    # Match sequences to libraries
    print("\nMatching sequences to libraries...")

    results = []
    library_matches = defaultdict(int)
    gene_counts = defaultdict(lambda: defaultdict(int))

    for seq in all_sequences:
        ctrl_count = ctrl_counts.get(seq, 0)
        loy_count = loy_counts.get(seq, 0)
        total_count = ctrl_count + loy_count

        record = {
            'sgRNA_sequence': seq,
            'D0_CTRL_count': ctrl_count,
            'D0_LOY_count': loy_count,
            'total_count': total_count,
            'library': 'Unknown',
            'gene': 'NA',
            'chromosome': 'NA',
            'start_pos': 'NA',
            'end_pos': 'NA',
            'strand': 'NA',
            'exon': 'NA'
        }

        # Check exact matches only
        if seq in avana_seqs:
            info = avana_seqs[seq]
            record.update({
                'library': 'Avana',
                'gene': info.get('gene', 'NA'),
                'chromosome': info.get('chromosome', 'NA'),
                'start_pos': info.get('start', 'NA'),
                'end_pos': info.get('end', 'NA'),
                'strand': info.get('strand', 'NA'),
                'exon': info.get('exon', 'NA')
            })
            library_matches['Avana'] += total_count
            gene_counts['Avana'][info.get('gene', 'NA')] += total_count

        elif seq in brunello_seqs:
            info = brunello_seqs[seq]
            record.update({
                'library': 'Brunello',
                'gene': info.get('gene', 'NA')
            })
            library_matches['Brunello'] += total_count
            gene_counts['Brunello'][info.get('gene', 'NA')] += total_count

        elif seq in tkov3_seqs:
            info = tkov3_seqs[seq]
            record.update({
                'library': 'TKOv3',
                'gene': info.get('gene', 'NA'),
                'chromosome': info.get('chromosome', 'NA'),
                'start_pos': info.get('start', 'NA'),
                'end_pos': info.get('end', 'NA'),
                'strand': info.get('strand', 'NA')
            })
            library_matches['TKOv3'] += total_count
            gene_counts['TKOv3'][info.get('gene', 'NA')] += total_count
        else:
            library_matches['Unknown'] += total_count

        results.append(record)

    # Create DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('total_count', ascending=False)

    # Save full results
    output_file = output_dir / "D0_unmapped_library_identification_with_coordinates.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nFull results saved to: {output_file}")

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)

    total_reads = sum(library_matches.values())
    print(f"\nTotal unmapped reads analyzed: {total_reads}")
    print(f"Unique sgRNA sequences: {len(all_sequences)}")

    print("\nLibrary Distribution (by read count):")
    print("-"*40)
    for lib in ['Avana', 'Brunello', 'TKOv3', 'Unknown']:
        count = library_matches[lib]
        pct = count / total_reads * 100 if total_reads > 0 else 0
        n_seqs = len(results_df[results_df['library'] == lib])
        print(f"{lib:10s}: {count:7d} reads ({pct:5.1f}%), {n_seqs:4d} unique sequences")

    # Top genes for each library
    for lib in ['Avana', 'Brunello', 'TKOv3']:
        if lib in gene_counts and len(gene_counts[lib]) > 0:
            print(f"\nTop 15 {lib} Genes (by read count):")
            print("-"*40)
            sorted_genes = sorted(gene_counts[lib].items(), key=lambda x: x[1], reverse=True)
            for i, (gene, count) in enumerate(sorted_genes[:15], 1):
                # Count unique sgRNAs for this gene
                n_sgrnas = len(results_df[(results_df['library'] == lib) & (results_df['gene'] == gene)])
                print(f"{i:2d}. {gene:20s}: {count:5d} reads ({n_sgrnas} unique sgRNAs)")

    # Show top unknown sequences
    unknown_df = results_df[results_df['library'] == 'Unknown'].head(20)
    if len(unknown_df) > 0:
        print("\nTop 20 Unknown Sequences (need BLAST):")
        print("-"*40)
        for i, row in enumerate(unknown_df.itertuples(), 1):
            print(f"{i:2d}. {row.sgRNA_sequence} - {row.total_count:3d} reads")

    # Create a summary with genomic distribution for Avana
    avana_df = results_df[results_df['library'] == 'Avana']
    if len(avana_df) > 0 and 'chromosome' in avana_df.columns:
        # Filter out NA chromosomes
        chrom_df = avana_df[avana_df['chromosome'] != 'NA']
        if len(chrom_df) > 0:
            print("\nAvana Contamination - Chromosomal Distribution:")
            print("-"*40)
            chrom_dist = chrom_df.groupby('chromosome')['total_count'].sum().sort_values(ascending=False)
            for chrom, count in chrom_dist.head(10).items():
                pct = count / chrom_df['total_count'].sum() * 100
                n_genes = len(chrom_df[chrom_df['chromosome'] == chrom]['gene'].unique())
                print(f"Chr {chrom:3s}: {count:5d} reads ({pct:5.1f}%), {n_genes} genes")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

    return results_df

if __name__ == "__main__":
    df = main()