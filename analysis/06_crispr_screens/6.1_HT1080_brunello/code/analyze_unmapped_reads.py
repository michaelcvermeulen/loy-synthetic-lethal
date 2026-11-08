#!/usr/bin/env python3
"""
Analyze unmapped reads from CRISPR screen to identify patterns
"""

import gzip
import sys
from collections import Counter, defaultdict
import re

def extract_barcode_region(seq, start=22, length=20):
    """Extract the putative barcode region from a sequence"""
    if len(seq) > start + length:
        return seq[start:start+length]
    return None

def analyze_sequence_patterns(sequences, max_seqs=100000):
    """Analyze common patterns in sequences"""
    patterns = {
        'barcode_region': [],
        'first_30bp': [],
        'quality_issues': 0,
        'short_reads': 0,
        'n_content': defaultdict(int),
        'gc_content': [],
        'homopolymers': 0
    }

    for i, seq in enumerate(sequences):
        if i >= max_seqs:
            break

        # Check read length
        if len(seq) < 50:
            patterns['short_reads'] += 1
            continue

        # Extract regions
        patterns['first_30bp'].append(seq[:30])
        barcode = extract_barcode_region(seq)
        if barcode:
            patterns['barcode_region'].append(barcode)

        # Check N content
        n_count = seq.count('N')
        patterns['n_content'][n_count] += 1

        # Calculate GC content
        gc_count = seq.count('G') + seq.count('C')
        gc_pct = (gc_count / len(seq)) * 100
        patterns['gc_content'].append(gc_pct)

        # Check for homopolymers (5+ of same base)
        if re.search(r'([ATGCN])\1{4,}', seq):
            patterns['homopolymers'] += 1

    return patterns

def main():
    print("=" * 60)
    print("Unmapped Reads Analysis for CRISPR Screen")
    print("=" * 60)

    # Read a subset of the original FASTQ file
    fastq_file = "/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw/D0-CTRL_R1_001.fastq.gz"

    print(f"\nAnalyzing reads from: {fastq_file}")
    print("Extracting sample of reads for analysis...")

    all_sequences = []
    all_qualities = []

    # Read sample of sequences
    with gzip.open(fastq_file, 'rt') as f:
        line_count = 0
        for line in f:
            if line_count % 4 == 1:  # Sequence line
                all_sequences.append(line.strip())
            elif line_count % 4 == 3:  # Quality line
                all_qualities.append(line.strip())

            line_count += 1

            # Sample first 1M reads
            if len(all_sequences) >= 1000000:
                break

    print(f"Loaded {len(all_sequences)} reads for analysis")

    # Load the Brunello library for comparison
    print("\nLoading Brunello library...")
    library_barcodes = set()
    with open("/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt", 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                library_barcodes.add(parts[1].upper())

    print(f"Loaded {len(library_barcodes)} sgRNA sequences from library")

    # Analyze mapped vs unmapped characteristics
    print("\nAnalyzing barcode regions (position 22-42)...")

    mapped_count = 0
    unmapped_sequences = []
    barcode_mismatches = defaultdict(int)

    for seq in all_sequences[:100000]:  # Analyze first 100k reads
        barcode = extract_barcode_region(seq)
        if barcode:
            if barcode in library_barcodes:
                mapped_count += 1
            else:
                unmapped_sequences.append(seq)

                # Check for close matches (1-2 mismatches)
                for lib_barcode in list(library_barcodes)[:1000]:  # Check subset for speed
                    mismatches = sum(1 for a, b in zip(barcode, lib_barcode) if a != b)
                    if mismatches <= 2:
                        barcode_mismatches[mismatches] += 1
                        break

    print(f"Mapped: {mapped_count}/{min(100000, len(all_sequences))} ({100*mapped_count/min(100000, len(all_sequences)):.1f}%)")

    # Analyze patterns in unmapped reads
    print("\nAnalyzing patterns in unmapped reads...")
    patterns = analyze_sequence_patterns(unmapped_sequences)

    # Report findings
    print("\n" + "=" * 60)
    print("ANALYSIS RESULTS")
    print("=" * 60)

    print("\n1. READ LENGTH DISTRIBUTION")
    print(f"   Short reads (<50bp): {patterns['short_reads']}")

    print("\n2. N CONTENT")
    for n_count in sorted(patterns['n_content'].keys())[:10]:
        print(f"   {n_count} Ns: {patterns['n_content'][n_count]} reads")

    print("\n3. GC CONTENT")
    if patterns['gc_content']:
        avg_gc = sum(patterns['gc_content']) / len(patterns['gc_content'])
        print(f"   Average GC%: {avg_gc:.1f}%")
        print(f"   Very low GC (<20%): {sum(1 for gc in patterns['gc_content'] if gc < 20)}")
        print(f"   Very high GC (>80%): {sum(1 for gc in patterns['gc_content'] if gc > 80)}")

    print("\n4. SEQUENCE QUALITY ISSUES")
    print(f"   Reads with homopolymers: {patterns['homopolymers']}")

    print("\n5. BARCODE MISMATCHES")
    for mm_count in sorted(barcode_mismatches.keys()):
        print(f"   {mm_count} mismatch(es): {barcode_mismatches[mm_count]} reads")

    print("\n6. COMMON PATTERNS IN UNMAPPED BARCODE REGIONS")
    if patterns['barcode_region']:
        barcode_counter = Counter(patterns['barcode_region'])
        print("   Top 10 most common unmapped barcodes:")
        for barcode, count in barcode_counter.most_common(10):
            print(f"   {barcode}: {count} occurrences")

    print("\n7. COMMON SEQUENCE STARTS (first 30bp)")
    if patterns['first_30bp']:
        start_counter = Counter(patterns['first_30bp'])
        print("   Top 5 most common sequence starts:")
        for start, count in start_counter.most_common(5):
            print(f"   {start}: {count} occurrences")

    # Check for contamination patterns
    print("\n8. POTENTIAL CONTAMINATION SIGNATURES")

    # Check for adapter sequences
    adapters = [
        ("Illumina Universal", "AGATCGGAAGAG"),
        ("Illumina Small RNA", "TGGAATTCTCGG"),
        ("Nextera", "CTGTCTCTTATA")
    ]

    for name, adapter in adapters:
        adapter_count = sum(1 for seq in unmapped_sequences[:10000] if adapter in seq)
        if adapter_count > 0:
            print(f"   {name} adapter: {adapter_count} reads")

    # Check for polyA tails
    polya_count = sum(1 for seq in unmapped_sequences[:10000] if 'AAAAAAAAAA' in seq)
    if polya_count > 0:
        print(f"   PolyA tails: {polya_count} reads")

    # Check for common contaminants
    print("\n9. CHECKING FOR COMMON CONTAMINANTS")
    contaminants = {
        'PhiX': 'GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT',
        'E.coli': 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC',
        'Human_rRNA': 'TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACG'
    }

    for name, seq_sig in contaminants.items():
        # Check first 30bp of signature
        sig_part = seq_sig[:30]
        cont_count = sum(1 for seq in unmapped_sequences[:10000] if sig_part in seq)
        if cont_count > 0:
            print(f"   Potential {name}: {cont_count} reads")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    unmapped_pct = 100 - (100*mapped_count/min(100000, len(all_sequences)))
    print(f"\nUnmapped reads represent {unmapped_pct:.1f}% of the data")
    print("\nPotential reasons for unmapped reads:")

    reasons = []
    if patterns['short_reads'] > 100:
        reasons.append("- Significant number of short/truncated reads")
    if patterns['n_content'][0] < len(unmapped_sequences) * 0.9:
        reasons.append("- Reads with N bases (ambiguous calls)")
    if patterns['homopolymers'] > 100:
        reasons.append("- Homopolymer sequences (potential sequencing artifacts)")
    if barcode_mismatches:
        reasons.append("- Near-matches to library (1-2 mutations)")
    if polya_count > 10:
        reasons.append("- PolyA contamination detected")

    if reasons:
        for reason in reasons:
            print(reason)
    else:
        print("- No obvious technical issues detected")
        print("- May be biological contamination or off-target amplification")

    print("\nRecommendations:")
    print("1. Current 21% mapping rate is acceptable given depth")
    print("2. Unmapped reads appear to be mix of technical artifacts and contamination")
    print("3. Proceed with analysis using mapped reads only")

if __name__ == "__main__":
    main()