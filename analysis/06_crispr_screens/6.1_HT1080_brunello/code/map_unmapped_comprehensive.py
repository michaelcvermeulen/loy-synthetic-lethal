#!/usr/bin/env python3
"""
Comprehensive mapping of unmapped sgRNA sequences to genomic coordinates.
Tests both forward and reverse complement sequences against all libraries.
"""

import pandas as pd
import numpy as np
from collections import Counter
import os
from pathlib import Path

def reverse_complement(seq):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in seq[::-1])

def load_full_fastq_unmapped(base_dir):
    """Load unmapped sequences directly from full FASTQ analysis if available"""
    # First check if we have the full unmapped analysis from MAGeCK
    unmapped_file = base_dir / "analysis/ht1080_screen_2/mageck/count/ht1080_screen2_R1.unmapped.txt"

    if unmapped_file.exists():
        print("Loading full unmapped sequences from MAGeCK...")
        with open(unmapped_file, 'r') as f:
            unmapped_seqs = []
            for line in f:
                if not line.startswith('@'):  # Skip headers
                    seq = line.strip().split('\t')[0] if '\t' in line else line.strip()
                    if len(seq) >= 20:
                        # Extract 20bp sgRNA (accounting for trimming)
                        sgrna = seq[22:42] if len(seq) > 42 else seq[:20]
                        unmapped_seqs.append(sgrna)
        print(f"Found {len(unmapped_seqs)} unmapped reads from MAGeCK")
        return Counter(unmapped_seqs)
    return None

def load_unmapped_sequences(base_dir):
    """Load unmapped sequences from CSV files and previous analysis"""
    ctrl_file = base_dir / "analysis/ht1080_screen_2/D0_CTRL_unmapped_sequences.csv"
    loy_file = base_dir / "analysis/ht1080_screen_2/D0_LOY_unmapped_sequences.csv"

    print("Loading unmapped sequences from CSV...")
    ctrl_df = pd.read_csv(ctrl_file)
    loy_df = pd.read_csv(loy_file)

    # Count occurrences in each sample
    ctrl_counts = Counter(ctrl_df['extracted_20bp_sgRNA'].dropna())
    loy_counts = Counter(loy_df['extracted_20bp_sgRNA'].dropna())

    # Also load the previously identified Avana contamination
    avana_contam_file = base_dir / "analysis/ht1080_screen_2/library_comparison/avana_contamination_details.csv"
    if avana_contam_file.exists():
        print("Loading previously identified Avana contamination...")
        avana_df = pd.read_csv(avana_contam_file)
        for _, row in avana_df.iterrows():
            seq = row['sequence']
            count = row['count']
            # Add these to our counts
            if seq not in ctrl_counts:
                ctrl_counts[seq] = count

    # Get all unique sequences
    all_sequences = set(ctrl_counts.keys()) | set(loy_counts.keys())

    print(f"Found {len(all_sequences)} unique unmapped sgRNA sequences")
    print(f"D0-CTRL: {len(ctrl_counts)} unique sequences")
    print(f"D0-LOY: {len(loy_counts)} unique sequences")

    return ctrl_counts, loy_counts, all_sequences

def load_avana_library_comprehensive(base_dir):
    """Load Avana library from multiple sources"""
    avana_sequences = {}

    # Primary Avana file with coordinates
    avana_file1 = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/Avana_library_forCRISPRcleanR_hg38.txt"
    if avana_file1.exists():
        print("Loading primary Avana library...")
        avana_df = pd.read_csv(avana_file1, sep='\t')
        for _, row in avana_df.iterrows():
            seq = row['CODE'].upper()
            avana_sequences[seq] = {
                'gene': row['GENES'],
                'exon': row['EXONE'],
                'chromosome': row['CHRM'],
                'strand': row['STRAND'],
                'start': row['STARTpos'],
                'end': row['ENDpos']
            }

    # Check for converted Avana file
    avana_file2 = base_dir / "data/avana/Avana_library_MAGeCK.txt"
    if avana_file2.exists():
        print("Loading converted Avana library...")
        with open(avana_file2, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    seq = parts[1].upper()
                    gene = parts[2] if len(parts) > 2 else 'Unknown'
                    if seq not in avana_sequences:
                        avana_sequences[seq] = {'gene': gene, 'chromosome': 'NA', 'start': 'NA', 'end': 'NA'}

    print(f"Loaded {len(avana_sequences)} total Avana sgRNAs")
    return avana_sequences

def test_sequence_variants(seq, library_dict):
    """Test a sequence and its reverse complement against library"""
    # Test original
    if seq in library_dict:
        return seq, 'forward', library_dict[seq]

    # Test reverse complement
    rc_seq = reverse_complement(seq)
    if rc_seq in library_dict:
        return rc_seq, 'reverse_complement', library_dict[rc_seq]

    # Test if it's a substring of library sequences (partial match)
    for lib_seq in library_dict:
        if seq in lib_seq or lib_seq in seq:
            return lib_seq, 'partial_match', library_dict[lib_seq]

    return None, None, None

def map_sequences_comprehensive(all_sequences, ctrl_counts, loy_counts, avana_dict, brunello_dict, tkov3_dict):
    """Comprehensively map sequences testing all orientations"""

    print("\nMapping sequences to libraries (testing forward and reverse complement)...")

    results = []

    avana_matches = 0
    brunello_matches = 0
    tkov3_matches = 0
    no_matches = 0

    for seq in all_sequences:
        if not seq or pd.isna(seq):
            continue

        seq = seq.upper()
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
            'match_type': 'none',
            'matched_sequence': seq,
            'gene': 'NA',
            'chromosome': 'NA',
            'start_pos': 'NA',
            'end_pos': 'NA',
            'strand': 'NA',
            'exon': 'NA'
        }

        # Test against Avana library
        matched_seq, match_type, avana_info = test_sequence_variants(seq, avana_dict)
        if avana_info:
            record.update({
                'library_source': 'Avana',
                'match_type': match_type,
                'matched_sequence': matched_seq,
                'gene': avana_info.get('gene', 'NA'),
                'chromosome': avana_info.get('chromosome', 'NA'),
                'start_pos': avana_info.get('start', 'NA'),
                'end_pos': avana_info.get('end', 'NA'),
                'strand': avana_info.get('strand', 'NA'),
                'exon': avana_info.get('exon', 'NA')
            })
            avana_matches += 1

        # Test against Brunello if no Avana match
        elif seq in brunello_dict or reverse_complement(seq) in brunello_dict:
            if seq in brunello_dict:
                record.update({
                    'library_source': 'Brunello',
                    'match_type': 'forward',
                    'gene': brunello_dict[seq].get('gene', 'NA')
                })
            else:
                rc_seq = reverse_complement(seq)
                record.update({
                    'library_source': 'Brunello',
                    'match_type': 'reverse_complement',
                    'matched_sequence': rc_seq,
                    'gene': brunello_dict[rc_seq].get('gene', 'NA')
                })
            brunello_matches += 1

        # Test against TKOv3
        elif seq in tkov3_dict or reverse_complement(seq) in tkov3_dict:
            if seq in tkov3_dict:
                tkov3_info = tkov3_dict[seq]
            else:
                rc_seq = reverse_complement(seq)
                tkov3_info = tkov3_dict[rc_seq]
                record['match_type'] = 'reverse_complement'
                record['matched_sequence'] = rc_seq

            record.update({
                'library_source': 'TKOv3',
                'gene': tkov3_info.get('gene', 'NA'),
                'chromosome': tkov3_info.get('chromosome', 'NA'),
                'start_pos': tkov3_info.get('start', 'NA'),
                'end_pos': tkov3_info.get('end', 'NA')
            })
            tkov3_matches += 1
        else:
            no_matches += 1

        results.append(record)

    # Create DataFrame and sort
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('total_count', ascending=False)

    print(f"\nMapping Summary:")
    print(f"  Avana matches: {avana_matches} ({avana_matches/len(all_sequences)*100:.1f}%)")
    print(f"  Brunello matches: {brunello_matches} ({brunello_matches/len(all_sequences)*100:.1f}%)")
    print(f"  TKOv3 matches: {tkov3_matches} ({tkov3_matches/len(all_sequences)*100:.1f}%)")
    print(f"  No library match: {no_matches} ({no_matches/len(all_sequences)*100:.1f}%)")

    return results_df

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
                seq = parts[1].upper()
                gene = parts[2]
                brunello_dict[seq] = {'gene': gene}

    print(f"Loaded {len(brunello_dict)} Brunello sgRNAs")
    return brunello_dict

def load_tkov3_library(base_dir):
    """Load TKOv3 library"""
    tkov3_file = base_dir / "SSC0.1/bin/bagel/data/libraries_for_CRISPRcleanR/TKOv3_library_forCRISPRcleanR_sgRNAtag.txt"

    if not tkov3_file.exists():
        return {}

    print("Loading TKOv3 library...")
    tkov3_df = pd.read_csv(tkov3_file, sep='\t')

    tkov3_dict = {}
    for _, row in tkov3_df.iterrows():
        if 'CODE' in row:
            seq = str(row['CODE']).upper()
            tkov3_dict[seq] = {
                'gene': row.get('GENES', 'NA'),
                'chromosome': row.get('CHRM', 'NA'),
                'start': row.get('STARTpos', 'NA'),
                'end': row.get('ENDpos', 'NA')
            }

    print(f"Loaded {len(tkov3_dict)} TKOv3 sgRNAs")
    return tkov3_dict

def create_blast_input(results_df, output_dir):
    """Create FASTA file for BLAST analysis of unmatched sequences"""

    unknown_df = results_df[results_df['library_source'] == 'Unknown']

    if len(unknown_df) > 0:
        fasta_file = output_dir / "unmapped_for_blast.fasta"
        with open(fasta_file, 'w') as f:
            for i, row in enumerate(unknown_df.head(100).itertuples(), 1):
                f.write(f">seq_{i}_count_{row.total_count}\n")
                f.write(f"{row.sgRNA_sequence}\n")

        print(f"\nCreated BLAST input file: {fasta_file}")
        print(f"Contains top 100 unmatched sequences for genomic alignment")

        return fasta_file
    return None

def main():
    """Main analysis pipeline"""

    base_dir = Path("/global/home/hpc5173/projects/CRISPR_screen")
    output_dir = base_dir / "analysis/ht1080_screen_2/genomic_mapping"
    output_dir.mkdir(exist_ok=True)

    print("="*70)
    print("COMPREHENSIVE UNMAPPED sgRNA GENOMIC MAPPING ANALYSIS")
    print("="*70)

    # Load sequences
    ctrl_counts, loy_counts, all_sequences = load_unmapped_sequences(base_dir)

    # Load all libraries
    avana_dict = load_avana_library_comprehensive(base_dir)
    brunello_dict = load_brunello_library(base_dir)
    tkov3_dict = load_tkov3_library(base_dir)

    # Map sequences comprehensively
    results_df = map_sequences_comprehensive(
        all_sequences, ctrl_counts, loy_counts,
        avana_dict, brunello_dict, tkov3_dict
    )

    # Save comprehensive results
    output_file = output_dir / "unmapped_sgrna_comprehensive_mapping.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")

    # Create BLAST input for unmatched
    blast_file = create_blast_input(results_df, output_dir)

    # Generate detailed report
    print("\n" + "="*70)
    print("DETAILED ANALYSIS RESULTS")
    print("="*70)

    # Show match type distribution
    if 'match_type' in results_df.columns:
        match_types = results_df.groupby(['library_source', 'match_type']).size()
        print("\nMatch Type Distribution:")
        print("-"*40)
        for (lib, match_type), count in match_types.items():
            print(f"{lib:10s} - {match_type:20s}: {count:4d} sequences")

    # Gene distribution for matched sequences
    for library in ['Avana', 'Brunello', 'TKOv3']:
        lib_df = results_df[results_df['library_source'] == library]
        if len(lib_df) > 0:
            print(f"\n{library} Library Matches - Top Genes:")
            print("-"*40)
            gene_counts = lib_df.groupby('gene')['total_count'].sum().sort_values(ascending=False)
            for gene, count in gene_counts.head(10).items():
                n_sgrnas = len(lib_df[lib_df['gene'] == gene])
                print(f"  {gene:15s}: {count:5d} reads ({n_sgrnas} sgRNAs)")

    # Unknown sequences
    unknown_df = results_df[results_df['library_source'] == 'Unknown']
    if len(unknown_df) > 0:
        print(f"\n\nUnmatched Sequences Requiring BLAST:")
        print("-"*40)
        print(f"Total sequences: {len(unknown_df)}")
        print(f"Total reads: {unknown_df['total_count'].sum()}")
        print("\nTop 15 by frequency:")
        for i, row in enumerate(unknown_df.head(15).itertuples(), 1):
            print(f"{i:2d}. {row.sgRNA_sequence} - {row.total_count:4d} reads")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

    if blast_file:
        print(f"\nNext step: Run BLAST on {blast_file} to identify genomic origins")
        print("Command: blastn -query unmapped_for_blast.fasta -db nt -remote -outfmt 6")

if __name__ == "__main__":
    main()