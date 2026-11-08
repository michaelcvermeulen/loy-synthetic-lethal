#!/usr/bin/env python3
"""
Map unmapped sgRNA sequences to genomic coordinates and identify source library.
This script aggregates all unmapped sequences, counts occurrences, and maps to known libraries.
"""

import pandas as pd
import numpy as np
from collections import Counter
import os
from pathlib import Path

def load_unmapped_sequences(base_dir):
    """Load unmapped sequences from CSV files"""
    ctrl_file = base_dir / "analysis/ht1080_screen_2/D0_CTRL_unmapped_sequences.csv"
    loy_file = base_dir / "analysis/ht1080_screen_2/D0_LOY_unmapped_sequences.csv"

    print("Loading unmapped sequences...")
    ctrl_df = pd.read_csv(ctrl_file)
    loy_df = pd.read_csv(loy_file)

    # Count occurrences in each sample
    ctrl_counts = Counter(ctrl_df['extracted_20bp_sgRNA'].dropna())
    loy_counts = Counter(loy_df['extracted_20bp_sgRNA'].dropna())

    # Get all unique sequences
    all_sequences = set(ctrl_counts.keys()) | set(loy_counts.keys())

    print(f"Found {len(all_sequences)} unique unmapped sgRNA sequences")
    print(f"D0-CTRL: {len(ctrl_counts)} unique sequences from {sum(ctrl_counts.values())} reads")
    print(f"D0-LOY: {len(loy_counts)} unique sequences from {sum(loy_counts.values())} reads")

    return ctrl_counts, loy_counts, all_sequences

def load_avana_library(base_dir):
    """Load Avana library with genomic coordinates"""
    avana_file = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/Avana_library_forCRISPRcleanR_hg38.txt"

    print("\nLoading Avana library with genomic coordinates...")
    avana_df = pd.read_csv(avana_file, sep='\t')

    # Create dictionary for fast lookup
    avana_dict = {}
    for _, row in avana_df.iterrows():
        seq = row['CODE']
        avana_dict[seq] = {
            'gene': row['GENES'],
            'exon': row['EXONE'],
            'chromosome': row['CHRM'],
            'strand': row['STRAND'],
            'start': row['STARTpos'],
            'end': row['ENDpos']
        }

    print(f"Loaded {len(avana_dict)} Avana sgRNAs with coordinates")
    return avana_dict

def load_brunello_library(base_dir):
    """Load Brunello library"""
    brunello_file = base_dir / "data/brunello/broadgpp-brunello-library-corrected.txt"

    print("Loading Brunello library...")
    brunello_dict = {}

    with open(brunello_file, 'r') as f:
        next(f)  # Skip header if present
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                seq = parts[1]
                gene = parts[2]
                brunello_dict[seq] = {'gene': gene}

    print(f"Loaded {len(brunello_dict)} Brunello sgRNAs")
    return brunello_dict

def load_tkov3_library(base_dir):
    """Load TKOv3 library if available"""
    tkov3_file = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/TKOv3_library_forCRISPRcleanR_sgRNAtag.txt"

    if not tkov3_file.exists():
        print("TKOv3 library file not found, skipping...")
        return {}

    print("Loading TKOv3 library...")
    tkov3_df = pd.read_csv(tkov3_file, sep='\t')

    tkov3_dict = {}
    for _, row in tkov3_df.iterrows():
        if 'CODE' in row and 'GENES' in row:
            seq = row['CODE']
            tkov3_dict[seq] = {
                'gene': row['GENES'],
                'chromosome': row.get('CHRM', 'NA'),
                'start': row.get('STARTpos', 'NA'),
                'end': row.get('ENDpos', 'NA')
            }

    print(f"Loaded {len(tkov3_dict)} TKOv3 sgRNAs")
    return tkov3_dict

def map_sequences_to_libraries(all_sequences, ctrl_counts, loy_counts, avana_dict, brunello_dict, tkov3_dict):
    """Map sequences to all available libraries and create comprehensive report"""

    print("\nMapping sequences to libraries...")

    results = []

    avana_matches = 0
    brunello_matches = 0
    tkov3_matches = 0
    no_matches = 0

    for seq in all_sequences:
        ctrl_count = ctrl_counts.get(seq, 0)
        loy_count = loy_counts.get(seq, 0)
        total_count = ctrl_count + loy_count

        # Initialize record
        record = {
            'sgRNA_sequence': seq,
            'D0_CTRL_count': ctrl_count,
            'D0_LOY_count': loy_count,
            'total_count': total_count,
            'library_source': 'Unknown',
            'gene': 'NA',
            'chromosome': 'NA',
            'start_pos': 'NA',
            'end_pos': 'NA',
            'strand': 'NA',
            'exon': 'NA'
        }

        # Check Avana library first (has coordinates)
        if seq in avana_dict:
            avana_info = avana_dict[seq]
            record.update({
                'library_source': 'Avana',
                'gene': avana_info['gene'],
                'chromosome': avana_info['chromosome'],
                'start_pos': avana_info['start'],
                'end_pos': avana_info['end'],
                'strand': avana_info['strand'],
                'exon': avana_info['exon']
            })
            avana_matches += 1

        # Check Brunello library
        elif seq in brunello_dict:
            record.update({
                'library_source': 'Brunello',
                'gene': brunello_dict[seq]['gene']
            })
            brunello_matches += 1

        # Check TKOv3 library
        elif seq in tkov3_dict:
            tkov3_info = tkov3_dict[seq]
            record.update({
                'library_source': 'TKOv3',
                'gene': tkov3_info['gene'],
                'chromosome': tkov3_info.get('chromosome', 'NA'),
                'start_pos': tkov3_info.get('start', 'NA'),
                'end_pos': tkov3_info.get('end', 'NA')
            })
            tkov3_matches += 1
        else:
            no_matches += 1

        results.append(record)

    # Create DataFrame and sort by total count
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('total_count', ascending=False)

    print(f"\nMapping Summary:")
    print(f"  Avana matches: {avana_matches} ({avana_matches/len(all_sequences)*100:.1f}%)")
    print(f"  Brunello matches: {brunello_matches} ({brunello_matches/len(all_sequences)*100:.1f}%)")
    print(f"  TKOv3 matches: {tkov3_matches} ({tkov3_matches/len(all_sequences)*100:.1f}%)")
    print(f"  No library match: {no_matches} ({no_matches/len(all_sequences)*100:.1f}%)")

    return results_df

def analyze_gene_distribution(results_df):
    """Analyze distribution of genes in contamination"""

    print("\n" + "="*60)
    print("GENE DISTRIBUTION ANALYSIS")
    print("="*60)

    # Filter for Avana matches only
    avana_df = results_df[results_df['library_source'] == 'Avana'].copy()

    if len(avana_df) > 0:
        # Group by gene and sum counts
        gene_counts = avana_df.groupby('gene')['total_count'].sum().sort_values(ascending=False)

        print("\nTop 20 Avana genes by total read count:")
        print("-"*50)
        for i, (gene, count) in enumerate(gene_counts.head(20).items(), 1):
            # Count how many different sgRNAs for this gene
            n_sgrnas = len(avana_df[avana_df['gene'] == gene])
            print(f"{i:2d}. {gene:15s} - {count:6d} reads ({n_sgrnas} sgRNAs)")

        # Chromosome distribution for Avana
        chrom_counts = avana_df.groupby('chromosome')['total_count'].sum().sort_values(ascending=False)
        print("\n\nChromosome distribution of Avana contamination:")
        print("-"*50)
        for chrom, count in chrom_counts.head(10).items():
            pct = count / avana_df['total_count'].sum() * 100
            print(f"Chr {chrom:3s}: {count:6d} reads ({pct:5.1f}%)")

    # Analyze unknown sequences
    unknown_df = results_df[results_df['library_source'] == 'Unknown'].copy()
    if len(unknown_df) > 0:
        print(f"\n\nUnknown sequences analysis:")
        print("-"*50)
        print(f"Total unique sequences: {len(unknown_df)}")
        print(f"Total read count: {unknown_df['total_count'].sum()}")
        print(f"Average reads per sequence: {unknown_df['total_count'].mean():.1f}")
        print(f"Median reads per sequence: {unknown_df['total_count'].median():.1f}")

        print("\nTop 10 unknown sequences by frequency:")
        for i, row in enumerate(unknown_df.head(10).itertuples(), 1):
            print(f"{i:2d}. {row.sgRNA_sequence} - {row.total_count:4d} reads")

def main():
    """Main analysis pipeline"""

    # Set up paths
    base_dir = Path("/global/home/hpc5173/projects/CRISPR_screen")
    output_dir = base_dir / "analysis/ht1080_screen_2/genomic_mapping"
    output_dir.mkdir(exist_ok=True)

    print("="*60)
    print("UNMAPPED sgRNA GENOMIC MAPPING ANALYSIS")
    print("="*60)

    # Load unmapped sequences
    ctrl_counts, loy_counts, all_sequences = load_unmapped_sequences(base_dir)

    # Load libraries
    avana_dict = load_avana_library(base_dir)
    brunello_dict = load_brunello_library(base_dir)
    tkov3_dict = load_tkov3_library(base_dir)

    # Map sequences to libraries
    results_df = map_sequences_to_libraries(
        all_sequences, ctrl_counts, loy_counts,
        avana_dict, brunello_dict, tkov3_dict
    )

    # Save results
    output_file = output_dir / "unmapped_sgrna_genomic_mapping.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")

    # Analyze gene distribution
    analyze_gene_distribution(results_df)

    # Create summary report
    summary_file = output_dir / "genomic_mapping_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("UNMAPPED sgRNA GENOMIC MAPPING SUMMARY\n")
        f.write("="*60 + "\n\n")

        f.write(f"Total unique unmapped sgRNAs: {len(all_sequences)}\n")
        f.write(f"Total unmapped reads analyzed: {sum(ctrl_counts.values()) + sum(loy_counts.values())}\n\n")

        library_stats = results_df.groupby('library_source').agg({
            'sgRNA_sequence': 'count',
            'total_count': 'sum'
        })

        f.write("Library Source Distribution:\n")
        f.write("-"*40 + "\n")
        for lib in ['Avana', 'Brunello', 'TKOv3', 'Unknown']:
            if lib in library_stats.index:
                n_seqs = library_stats.loc[lib, 'sgRNA_sequence']
                n_reads = library_stats.loc[lib, 'total_count']
                pct_seqs = n_seqs / len(all_sequences) * 100
                pct_reads = n_reads / results_df['total_count'].sum() * 100
                f.write(f"{lib:10s}: {n_seqs:4d} sequences ({pct_seqs:5.1f}%), {n_reads:6d} reads ({pct_reads:5.1f}%)\n")

        # Add gene summary for Avana
        avana_df = results_df[results_df['library_source'] == 'Avana']
        if len(avana_df) > 0:
            f.write("\n\nTop Avana Contaminating Genes:\n")
            f.write("-"*40 + "\n")
            gene_counts = avana_df.groupby('gene')['total_count'].sum().sort_values(ascending=False)
            for gene, count in gene_counts.head(20).items():
                f.write(f"{gene:15s}: {count:6d} reads\n")

    print(f"\nSummary saved to: {summary_file}")
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()