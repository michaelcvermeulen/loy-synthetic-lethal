#!/usr/bin/env python3
"""
Align unknown sgRNA sequences to human genome (hg38) using Biopython.
Identifies genomic origins of contaminating sequences.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import time
import re

# Biopython imports
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

# Set up Entrez
Entrez.email = "crispr.screen.analysis@example.com"

def load_unknown_sequences(mapping_file):
    """Load unknown sequences from previous mapping analysis"""

    print("Loading unknown sequences...")
    df = pd.read_csv(mapping_file)

    # Filter for unknown sequences only
    unknown_df = df[df['library'] == 'Unknown'].copy()
    unknown_df = unknown_df.sort_values('total_count', ascending=False)

    print(f"Found {len(unknown_df)} unknown sequences")
    print(f"Total reads for unknown sequences: {unknown_df['total_count'].sum()}")

    return unknown_df

def create_sequence_records(unknown_df, max_sequences=100):
    """Create BioPython SeqRecord objects for BLAST"""

    records = []

    # Take top sequences by frequency
    for idx, row in enumerate(unknown_df.head(max_sequences).itertuples(), 1):
        seq = Seq(row.sgRNA_sequence)

        # Create record with metadata
        record = SeqRecord(
            seq,
            id=f"unknown_{idx}",
            description=f"count={row.total_count} ctrl={row.D0_CTRL_count} loy={row.D0_LOY_count}"
        )
        records.append(record)

        # Also add reverse complement
        rc_record = SeqRecord(
            seq.reverse_complement(),
            id=f"unknown_{idx}_rc",
            description=f"reverse_complement count={row.total_count}"
        )
        records.append(rc_record)

    print(f"Created {len(records)} sequence records (including reverse complements)")
    return records

def save_fasta_for_blast(records, output_file):
    """Save sequences in FASTA format for BLAST"""

    SeqIO.write(records, output_file, "fasta")
    print(f"Saved {len(records)} sequences to {output_file}")
    return output_file

def check_pam_sequences(sequence, pam="GG"):
    """Check if sequence has adjacent PAM site (NGG for SpCas9)"""

    # Check if PAM is immediately downstream
    # For CRISPR, PAM should be at positions 21-23 (NGG)
    # But since these are 20bp, check if they end with GG

    pam_present = False
    pam_position = None

    if sequence.endswith(pam):
        pam_present = True
        pam_position = "3_prime_end"

    return pam_present, pam_position

def analyze_sequence_properties(unknown_df):
    """Analyze properties of unknown sequences"""

    print("\nAnalyzing sequence properties...")

    properties = []

    for _, row in unknown_df.iterrows():
        seq = row['sgRNA_sequence']

        # Calculate properties
        gc = gc_fraction(seq)
        length = len(seq)

        # Check for PAM
        pam_present, pam_pos = check_pam_sequences(seq)

        # Check for homopolymers
        has_polyA = 'AAAA' in seq
        has_polyT = 'TTTT' in seq
        has_polyG = 'GGGG' in seq
        has_polyC = 'CCCC' in seq
        has_homopolymer = any([has_polyA, has_polyT, has_polyG, has_polyC])

        # Check for low complexity
        unique_kmers = len(set([seq[i:i+3] for i in range(len(seq)-2)]))
        low_complexity = unique_kmers < 10

        properties.append({
            'sequence': seq,
            'count': row['total_count'],
            'gc_content': gc,
            'length': length,
            'has_pam': pam_present,
            'has_homopolymer': has_homopolymer,
            'low_complexity': low_complexity,
            'unique_3mers': unique_kmers
        })

    props_df = pd.DataFrame(properties)

    # Summary statistics
    print(f"Average GC content: {props_df['gc_content'].mean():.2%}")
    print(f"Sequences with PAM: {props_df['has_pam'].sum()} ({props_df['has_pam'].mean():.1%})")
    print(f"Sequences with homopolymers: {props_df['has_homopolymer'].sum()} ({props_df['has_homopolymer'].mean():.1%})")
    print(f"Low complexity sequences: {props_df['low_complexity'].sum()} ({props_df['low_complexity'].mean():.1%})")

    return props_df

def run_local_blast(fasta_file, db_path=None):
    """Run local BLAST if database is available"""

    import subprocess
    import os

    # Check if BLAST is available
    try:
        result = subprocess.run(['which', 'blastn'], capture_output=True, text=True)
        if result.returncode != 0:
            print("BLAST not found in PATH, skipping local alignment")
            return None
    except:
        return None

    # Run BLAST
    output_file = fasta_file.replace('.fasta', '_blast_results.xml')

    cmd = [
        'blastn',
        '-query', fasta_file,
        '-db', 'nt',  # or specify local database
        '-out', output_file,
        '-outfmt', '5',  # XML output
        '-max_target_seqs', '10',
        '-word_size', '7',  # Shorter word size for short sequences
        '-evalue', '1000',  # Relaxed E-value for short sequences
        '-remote'  # Use remote BLAST
    ]

    print(f"Running BLAST (this may take several minutes)...")
    print(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode == 0:
            print(f"BLAST completed successfully. Results saved to {output_file}")
            return output_file
        else:
            print(f"BLAST failed: {result.stderr}")
            return None
    except subprocess.TimeoutExpired:
        print("BLAST timed out after 10 minutes")
        return None
    except Exception as e:
        print(f"Error running BLAST: {e}")
        return None

def run_ncbi_blast_batch(sequences, batch_size=10):
    """Run NCBI BLAST in batches to avoid timeouts"""

    results = []

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i+batch_size]
        print(f"\nProcessing batch {i//batch_size + 1}/{(len(sequences)-1)//batch_size + 1}")

        # Create query sequence
        query = '\n'.join([f">{s.id}\n{str(s.seq)}" for s in batch])

        try:
            print("Submitting to NCBI BLAST...")
            result_handle = NCBIWWW.qblast(
                "blastn",
                "nt",
                query,
                entrez_query="Homo sapiens[Organism]",
                megablast=True,
                hitlist_size=5,
                expect=1000,  # Relaxed E-value for short sequences
                word_size=7
            )

            print("Parsing BLAST results...")
            blast_records = NCBIXML.parse(result_handle)

            for record in blast_records:
                results.append(record)

            # Be nice to NCBI
            time.sleep(3)

        except Exception as e:
            print(f"Error in BLAST batch: {e}")
            continue

    return results

def parse_blast_results(blast_file):
    """Parse BLAST XML results"""

    print(f"\nParsing BLAST results from {blast_file}...")

    alignments = []

    with open(blast_file, 'r') as handle:
        blast_records = NCBIXML.parse(handle)

        for record in blast_records:
            query_id = record.query
            query_length = record.query_length

            if len(record.alignments) == 0:
                alignments.append({
                    'query_id': query_id,
                    'query_length': query_length,
                    'hits': 0,
                    'best_hit': 'No hits',
                    'best_evalue': None,
                    'best_identity': None,
                    'chromosome': None,
                    'position': None
                })
            else:
                # Process top hit
                top_alignment = record.alignments[0]
                top_hsp = top_alignment.hsps[0]

                # Extract chromosome from title if present
                chromosome = extract_chromosome(top_alignment.title)

                # Calculate percent identity
                identity = (top_hsp.identities / top_hsp.align_length) * 100

                alignments.append({
                    'query_id': query_id,
                    'query_length': query_length,
                    'hits': len(record.alignments),
                    'best_hit': top_alignment.title[:100],
                    'best_evalue': top_hsp.expect,
                    'best_identity': identity,
                    'best_score': top_hsp.score,
                    'chromosome': chromosome,
                    'subject_start': top_hsp.sbjct_start,
                    'subject_end': top_hsp.sbjct_end,
                    'query_seq': str(top_hsp.query),
                    'match_seq': str(top_hsp.match),
                    'subject_seq': str(top_hsp.sbjct)
                })

    results_df = pd.DataFrame(alignments)

    # Summary
    print(f"Total queries: {len(results_df)}")
    print(f"Queries with hits: {(results_df['hits'] > 0).sum()}")
    print(f"Average best identity: {results_df['best_identity'].mean():.1f}%")

    return results_df

def extract_chromosome(title):
    """Extract chromosome from BLAST hit title"""

    # Common patterns in NCBI titles
    patterns = [
        r'chromosome (\d+|[XY])',
        r'chr(\d+|[XY])',
        r'Homo sapiens chromosome (\d+|[XY])',
    ]

    for pattern in patterns:
        match = re.search(pattern, title, re.IGNORECASE)
        if match:
            return match.group(1)

    return None

def categorize_genomic_hits(blast_results_df):
    """Categorize BLAST hits by genomic feature type"""

    print("\nCategorizing genomic hits...")

    categories = defaultdict(list)

    for _, row in blast_results_df.iterrows():
        if row['hits'] == 0:
            categories['no_hit'].append(row['query_id'])
        elif row['best_identity'] > 95:
            categories['high_identity'].append(row['query_id'])
        elif row['best_identity'] > 85:
            categories['moderate_identity'].append(row['query_id'])
        else:
            categories['low_identity'].append(row['query_id'])

    # Chromosome distribution
    if 'chromosome' in blast_results_df.columns:
        chrom_dist = blast_results_df['chromosome'].value_counts()
        print("\nChromosome distribution of hits:")
        for chrom, count in chrom_dist.head(10).items():
            if chrom:
                print(f"  Chr {chrom}: {count} hits")

    print("\nHit categories:")
    for category, queries in categories.items():
        print(f"  {category}: {len(queries)} sequences")

    return categories

def generate_report(unknown_df, props_df, blast_results_df, output_dir):
    """Generate comprehensive report of genomic alignment"""

    report_file = output_dir / "unknown_sequences_hg38_alignment_report.md"

    with open(report_file, 'w') as f:
        f.write("# Unknown sgRNA Sequences - hg38 Genome Alignment Report\n\n")

        f.write("## Executive Summary\n\n")
        f.write(f"- Analyzed {len(unknown_df)} unknown sgRNA sequences\n")
        f.write(f"- Total reads represented: {unknown_df['total_count'].sum()}\n")

        if blast_results_df is not None and len(blast_results_df) > 0:
            hits = (blast_results_df['hits'] > 0).sum()
            pct_hits = hits / len(blast_results_df) * 100
            f.write(f"- Sequences with genomic hits: {hits}/{len(blast_results_df)} ({pct_hits:.1f}%)\n")

            if 'best_identity' in blast_results_df.columns:
                avg_identity = blast_results_df[blast_results_df['hits'] > 0]['best_identity'].mean()
                f.write(f"- Average sequence identity: {avg_identity:.1f}%\n")

        f.write("\n## Sequence Properties\n\n")
        f.write(f"- Average GC content: {props_df['gc_content'].mean():.2%}\n")
        f.write(f"- Sequences with PAM (NGG): {props_df['has_pam'].sum()} ({props_df['has_pam'].mean():.1%})\n")
        f.write(f"- Sequences with homopolymers: {props_df['has_homopolymer'].sum()}\n")
        f.write(f"- Low complexity sequences: {props_df['low_complexity'].sum()}\n")

        f.write("\n## Top Unknown Sequences\n\n")
        f.write("| Rank | Sequence | Count | GC% | PAM | Identity% | Hit |\n")
        f.write("|------|----------|-------|-----|-----|-----------|-----|\n")

        for i, row in enumerate(unknown_df.head(20).itertuples(), 1):
            seq = row.sgRNA_sequence
            gc = gc_fraction(seq)
            pam, _ = check_pam_sequences(seq)

            # Get BLAST result if available
            if blast_results_df is not None and len(blast_results_df) > 0:
                blast_info = blast_results_df[blast_results_df['query_id'].str.contains(f"unknown_{i}")].head(1)
                if len(blast_info) > 0:
                    identity = blast_info.iloc[0]['best_identity']
                    hit = "Yes" if blast_info.iloc[0]['hits'] > 0 else "No"
                else:
                    identity = None
                    hit = "N/A"
            else:
                identity = None
                hit = "N/A"

            f.write(f"| {i} | {seq} | {row.total_count} | {gc:.1%} | {'Yes' if pam else 'No'} | ")
            if identity is not None:
                f.write(f"{identity:.1f} | {hit} |\n")
            else:
                f.write(f"N/A | {hit} |\n")

        if blast_results_df is not None and len(blast_results_df) > 0:
            f.write("\n## Genomic Alignment Results\n\n")

            # Chromosome distribution
            if 'chromosome' in blast_results_df.columns:
                chrom_counts = blast_results_df['chromosome'].value_counts()
                f.write("### Chromosome Distribution\n\n")
                f.write("| Chromosome | Hits |\n")
                f.write("|------------|------|\n")
                for chrom, count in chrom_counts.head(10).items():
                    if chrom:
                        f.write(f"| {chrom} | {count} |\n")

            f.write("\n### Example High-Identity Hits\n\n")
            high_identity = blast_results_df[blast_results_df['best_identity'] > 90].head(5)
            for _, hit in high_identity.iterrows():
                f.write(f"**{hit['query_id']}**\n")
                f.write(f"- Best hit: {hit['best_hit'][:80]}...\n")
                f.write(f"- Identity: {hit['best_identity']:.1f}%\n")
                f.write(f"- E-value: {hit['best_evalue']:.2e}\n\n")

        f.write("\n## Conclusions\n\n")
        f.write("Based on the analysis:\n\n")

        if blast_results_df is not None and (blast_results_df['hits'] > 0).sum() > len(blast_results_df) * 0.5:
            f.write("1. **Genomic DNA contamination likely** - Majority of unknown sequences align to human genome\n")
        else:
            f.write("1. **Technical artifacts or degraded library** - Most sequences don't align well to genome\n")

        if props_df['has_pam'].mean() > 0.2:
            f.write("2. **Potential off-target sites** - Significant fraction have PAM sequences\n")

        if props_df['low_complexity'].mean() > 0.3:
            f.write("3. **PCR artifacts possible** - Many low-complexity sequences detected\n")

        f.write("\n---\n")
        f.write("Report generated using Biopython-enhanced analysis pipeline\n")

    print(f"\nReport saved to: {report_file}")
    return report_file

def main():
    """Main analysis pipeline"""

    base_dir = Path("/global/home/hpc5173/projects/CRISPR_screen")
    output_dir = base_dir / "analysis/ht1080_screen_2/genomic_mapping"

    print("="*70)
    print("UNKNOWN sgRNA SEQUENCES - hg38 GENOME ALIGNMENT")
    print("="*70)

    # Load unknown sequences from previous analysis
    mapping_file = output_dir / "D0_unmapped_library_identification_with_coordinates.csv"
    unknown_df = load_unknown_sequences(mapping_file)

    # Analyze sequence properties
    props_df = analyze_sequence_properties(unknown_df)

    # Create sequence records for BLAST
    records = create_sequence_records(unknown_df, max_sequences=50)  # Start with top 50

    # Save as FASTA
    fasta_file = output_dir / "unknown_sequences_for_blast.fasta"
    save_fasta_for_blast(records, fasta_file)

    # Try local BLAST first
    blast_results_file = run_local_blast(str(fasta_file))

    if blast_results_file and Path(blast_results_file).exists():
        # Parse BLAST results
        blast_results_df = parse_blast_results(blast_results_file)

        # Categorize hits
        categories = categorize_genomic_hits(blast_results_df)

        # Save detailed results
        blast_results_df.to_csv(output_dir / "blast_results_detailed.csv", index=False)
    else:
        print("\nLocal/remote BLAST not available or failed")
        print("Consider running BLAST manually or using NCBI web interface")
        blast_results_df = None

    # Generate comprehensive report
    report_file = generate_report(unknown_df, props_df, blast_results_df, output_dir)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nKey outputs:")
    print(f"1. FASTA file: {fasta_file}")
    print(f"2. Report: {report_file}")
    if blast_results_df is not None:
        print(f"3. BLAST results: {output_dir}/blast_results_detailed.csv")

    return unknown_df, props_df, blast_results_df

if __name__ == "__main__":
    unknown_df, props_df, blast_results_df = main()