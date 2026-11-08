"""Assemble a standalone report from completed stage outputs; no heavy analysis."""
import base64
from datetime import datetime
from html import escape
import json
from pathlib import Path
import pandas as pd
from locus_plots import create_loci

ROOT = Path(__file__).resolve().parents[1]


def table(frame):
    return frame.to_html(index=False, border=0, float_format=lambda x:f'{x:.4g}', escape=True)


def figure(path, caption):
    full = ROOT/path
    if not full.exists():
        return ''
    encoded = base64.b64encode(full.read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{encoded}" alt="{escape(caption)}"><figcaption>{escape(caption)}</figcaption></figure>'


if __name__ == '__main__':
    samples = pd.read_csv(ROOT/'config/samples.tsv',sep='\t')
    qc = pd.read_csv(ROOT/'qc/qc_metrics_master.csv')
    consensus = pd.read_csv(ROOT/'peaks/consensus_membership.tsv',sep='\t')
    summaries = pd.read_csv(ROOT/'differential/analysis_summary.csv')
    primary = pd.read_csv(ROOT/'differential/primary_TMM_all.csv')
    annotated = pd.read_csv(ROOT/'annotation/annotated_differential_peaks.csv')
    create_loci(ROOT, annotated, samples)
    primary_summary = summaries.loc[summaries.arm == 'primary_TMM'].iloc[0]
    go_summaries = []
    for direction in ['gained', 'lost']:
        path = ROOT/f'annotation/GO_BP_promoter_{direction}.csv'
        if path.exists():
            go = pd.read_csv(path)
            go_summaries.append({'direction':direction, 'reported_terms':len(go),
                'terms_BH_FDR_below_005':int((go['p.adjust'] < .05).sum()),
                'minimum_adjusted_P':go['p.adjust'].min()})
    pieces = ['<!doctype html><html><head><meta charset="utf-8"><title>HT1080 ATACseq reprocessing</title>',
      '<style>body{font:16px/1.55 system-ui,sans-serif;max-width:1150px;margin:40px auto;padding:0 24px;color:#18242b}'
      'h1,h2{line-height:1.2}h2{margin-top:40px}table{border-collapse:collapse;width:100%;font-size:13px;display:block;overflow:auto}'
      'td,th{padding:8px 12px;border-bottom:1px solid #d9e1e5;text-align:left}th{background:#edf3f5}'
      'figure{margin:28px 0}img{max-width:100%;height:auto}figcaption,.muted{color:#52646d;font-size:14px}'
      '.note{background:#f0f5f7;border-left:4px solid #46788e;padding:14px 18px}code{background:#edf3f5;padding:2px 4px}</style></head><body>',
      '<h1>HT1080 ATACseq: LOY versus CTRL</h1>',
      f'<p class="muted">Run: SEPT_08_2026_REPROCESS. Report assembled {datetime.now().isoformat(timespec="seconds")}.</p>',
      '<p>Fresh peaks and accessibility analysis from six existing hg38 BAMs. The primary design treats three CTRL and three LOY samples as independent replicates, per the user’s instruction. Positive log2 fold changes indicate greater accessibility in LOY.</p>',
      f'<p><b>{len(consensus):,} consensus regions</b> were constructed. The primary model tested <b>{int(primary_summary.tested):,}</b> regions after excluding chrY/PAR and low-count regions, with <b>{int(primary_summary.gained):,} gained</b> and <b>{int(primary_summary.lost):,} lost</b> at FDR &lt;0.05.</p>',
      '<h2>Samples and quality</h2>',table(qc[['sample','accepted_fragments','tss_enrichment','nfr_fraction','n_peaks','frip_own','frip_consensus','unmatched_primary_records']]),
      '<p>TSS enrichment is computed from +4/−5 corrected insertion coordinates, using the mean in the central 100 bp divided by the outer 100 bp flanks of a ±2 kb window. The pinned knownGene TSS list excludes chrY, duplicate TSS positions, boundary-truncated windows, and blacklist-overlapping windows. FRiP uses unique eligible fragments overlapping each peak union. Original mitochondrial contamination and raw alignment rate cannot be recovered from the already filtered BAMs.</p>',
      figure('qc/tss_profiles.png','TSS insertion enrichment, computed consistently across all samples.'),
      figure('qc/fragment_lengths.png','Fragment length distributions; values are fractions of retained fragments.'),
      figure('qc/pca.png','PCA of normalized logCPM over the tested peak set.'),
      figure('qc/mds.png','MDS of normalized peak counts.'),
      figure('qc/correlation_heatmap.png','Spearman correlations over the common tested peak set.'),
      figure('qc/chrY_signal.png','chrY signal relative to eligible non-Y fragments. This is an ATAC signal check, not an estimate of the fraction of Y-negative cells.'),
      '<h2>Differential accessibility and sensitivity</h2>',table(summaries),
      '<p>Primary analysis: edgeR quasi-likelihood with peak-based TMM and design <code>~ condition</code>. Sensitivity arms use peak-free 10 kb background-bin normalization or passage adjustment. All arms test the same abundance-filtered peak matrix. The TREAT column tests effects beyond |log2 fold change|=0.5 and controls its own BH FDR.</p>',
      '<p>The passage-adjusted model produces many more significant regions while retaining very similar effect estimates. Statistical significance is therefore sensitive to the treatment of passage-associated variation. The independent-replicate design remains the primary analysis, as specified by the user.</p>',
      '<div class="note">Normalization differences are a limitation of interpretation. Peak counts, FRiP, and TSS enrichment do not establish an absolute genome-wide accessibility shift. Background-bin normalization is a sensitivity analysis; copy-number stability and external calibration have not yet been established.</div>',
      figure('differential/primary_volcano.png','Primary edgeR results; positive effects denote LOY gains.'),
      figure('qc/gc_effects.png','GC dependence of estimated accessibility changes under the primary and background normalization arms. This diagnostic alone cannot separate biological from technical GC effects.'),
      '<h2>Annotated regions</h2>',
      '<p>Promoter overlap uses ±2 kb around the pinned TSS annotation. All promoter-to-gene links are retained in the tables. Nearest-gene assignments for distal regions are hypotheses. The FDR values below remain peak-level statistics.</p>',
      table(annotated.sort_values('FDR')[['peak_id','annotation','nearest_symbol','logFC','FDR','direction']].head(25)),
      figure('annotation/example_loci.png','Two top primary gained and two top primary lost loci, separated by at least 10 kb where available. All replicates are shown. Coverage uses non-Y fragment depth scaling, not the model TMM factors; shaded intervals are tested peaks. Nearest genes are positional assignments.'),
      '<p>Promoter gene summaries contain median peak effects and gained/lost peak counts. Minimum peak FDR is descriptive and is not a gene-level FDR.</p>',
      table(pd.DataFrame(go_summaries)),
      '<p>GO figures rank terms by adjusted P; ranked terms are not necessarily significant. Small foreground gene counts and related, overlapping GO categories warrant cautious interpretation.</p>',
      figure('annotation/GO_BP_promoter_gained.png','Promoter-associated GO biological processes for gained peaks. The universe is genes with tested promoter peaks; remaining peak-number bias limits interpretation.'),
      figure('annotation/GO_BP_promoter_lost.png','Promoter-associated GO biological processes for lost peaks, with the same tested-gene universe.')]
    support_file = ROOT/'differential/primary_TMM_csaw_support.csv'
    if support_file.exists():
        support = pd.read_csv(support_file)
        number = ((support.FDR < 0.05) & support.csaw_significant_same_direction).sum()
        pieces += ['<h2>Window-based corroboration</h2>',
          f'<p>{number:,} primary significant peaks overlap a significant csaw region with the same representative effect direction. csaw uses 150 bp windows, region-level FDR after merging, and the background normalization sensitivity arm. The shared reads mean this is method robustness, not independent experimental validation.</p>']
    depth_file = ROOT/'qc/matched_depth_peak_counts.tsv'
    if depth_file.exists():
        pieces += ['<h2>Matched-depth peak discovery</h2>',table(pd.read_csv(depth_file,sep='\t')),
          '<p>Each sample was sampled without replacement to exactly the same number of eligible fragments using a fixed seed, then peaks were recalled with matching MACS3 settings. The primary analysis retains full-depth data.</p>']
    motif_file = ROOT/'motifs/motif_enrichment.csv'
    if motif_file.exists():
        motifs = pd.read_csv(motif_file)
        pieces += ['<h2>Original motifmatchr/JASPAR2022 enrichment screen</h2>',
          '<p>The original screen used motifmatchr with human JASPAR2022 CORE motifs and custom stratified enrichment tests in gained/lost peaks against nonsignificant accessible controls matched on GC content and accessibility. Tests were adjusted across all motifs and both directions. Matching losses are recorded. The HOMER/JASPAR2024 extension below uses the additional |log2FC| &gt;1 peak filter and a different enrichment procedure, so these result counts are not directly comparable. Neither screen demonstrates TF binding or activity.</p>',
          table(pd.read_csv(ROOT/'motifs/matching_accounting.csv')),
          table(motifs.head(20)),figure('motifs/gained_motifs.png','Motif enrichment in gained regions.'),
          figure('motifs/lost_motifs.png','Motif enrichment in lost regions.')]
    tf_dir = ROOT/'motifs/TF_extension'
    if (tf_dir/'chromVAR/analysis_summary.csv').exists():
        cv_summary = pd.read_csv(tf_dir/'chromVAR/analysis_summary.csv')
        cv = pd.read_csv(tf_dir/'chromVAR/motifs_independent_limma.csv')
        pieces += ['<h2 id="tf-extension">TF analysis extension: per-sample motif accessibility</h2>',
          '<p>chromVAR uses the six separate integer fragment-count columns and JASPAR2024 CORE vertebrate motifs. The peak universe excludes chrY/PAR and four sequences containing N bases. Background sampling matches GC and average accessibility. No DA threshold selects the peaks entering chromVAR.</p>',
          '<p>Technical batch is not recorded in the current ATAC metadata. Passage labels 11/12/21 are retained as a separate sensitivity term; they have not been relabeled as technical batch. Independent-replicate results remain the current primary comparison.</p>',
          table(cv_summary),
          '<p>Limma compares sample-level chromVAR z-scores with BH correction within each analysis arm. The effect is a LOY-minus-CTRL z-score difference, not log2 fold change. Welch t-tests are provided separately. The much larger number of passage-adjusted calls makes design sensitivity consequential.</p>',
          figure('motifs/TF_extension/chromVAR/motifs_independent_heatmap.png','Top-ranked individual motifs. Values are chromVAR deviation z-scores without row standardization. Ranked rows need not pass FDR <0.05.'),
          figure('motifs/TF_extension/chromVAR/archetypes_independent_heatmap.png','Vierstra-mapped motif groups use unions of member motif-positive peaks. Ambiguous or unmapped motifs remain individual entries labeled unmapped_; the archetype analysis therefore includes both collapsed groups and unresolved individual motifs.'),
          table(cv[['id','motif_name','LOY_minus_CTRL_zscore','adj.P.Val','archetype','mapping_status']].head(15)),
          '<p>Exact JASPAR-version mappings are preferred; unique accession matches to other versions are labeled as such. The full mapping audit and member identities are saved. Family summaries do not establish which paralog binds.</p>',
          figure('motifs/TF_extension/concordance/chromVAR_TF_RNA.png','TF expression from the selected HT1080 batch + treatment RNA result compared with motif accessibility. Exact and case-insensitive candidate symbol links are distinguished in the table; these are group effects, not matched-sample correlations.'),
          '<p>The leading independent-model motif changes include GATA- and FOX-related motifs with lower accessibility in LOY. Transcript direction is not a required agreement criterion; for example, GATA4 RNA increases while its chromVAR motif score decreases. Motif accessibility is not a direct assay of TF protein activity or binding.</p>',
          table(pd.read_csv(tf_dir/'concordance/sex_chromosome_TF_coverage.csv')),
          '<p>Selected sex-chromosome TF coverage: ZFY RNA decreases in the selected batch-adjusted RNA result, but this JASPAR set has no motif named ZFY. Zfx is a case-insensitive candidate label match to ZFX; a motif label does not establish paralog-specific binding. SRY has a motif model but is absent from the selected RNA result table. Missing labels or RNA rows do not establish absence of protein expression or occupancy.</p>',
          '<h3>Requested DA effect-size subset</h3>',
          '<p>The supplementary DA plots highlight edgeR FDR <0.05 and |log2FC| >1. This effect-size filter does not replace a minimum-effect hypothesis test. Known-motif enrichment and footprinting stages are tracked separately below.</p>',
          figure('motifs/TF_extension/HOMER/primary_TMM/DA_volcano_MA.png','Independent-replicate DA volcano and MA plots at the requested cutoffs.'),
          '<p>Methods, score matrices, all tests, motif mappings, and TF expression links are under <code>motifs/TF_extension/</code>; execution choices are documented in <code>TF_EXTENSION_PLAN.md</code>.</p>',
          '<p>Published methods: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC2898526/">HOMER (Heinz et al., 2010)</a>; <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5623146/">chromVAR (Schep et al., 2017)</a>; <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7449963/">TOBIAS (Bentsen et al., 2020)</a>. Published tools support the analysis framework; the conclusions also depend on this experiment\'s replication, matching, design, and QC.</p>']
        homer_summary = tf_dir/'HOMER/known_motif_summary.csv'
        if homer_summary.exists():
            pieces += ['<h3>HOMER known-motif enrichment</h3>',
              '<p>HOMER 5.1 tests the same 879 JASPAR2024 motifs with exact peak widths against GC/accessibility-matched non-DA consensus regions. The normalized hg38 reference and every extracted peak sequence pass independent checks. P-values use the more precise log-P output; BH correction spans both directions within each model arm. Fold enrichment uses HOMER motif fractions after its background weighting.</p>',
              table(pd.read_csv(tf_dir/'HOMER/analysis_summary.csv')),table(pd.read_csv(homer_summary).head(24)),
              '<p>Six motif models pass FDR &lt;0.05 in both the primary lost-peak HOMER analysis and independent-replicate chromVAR comparison: GATA6, TRPS1, Gata3, FOXD2, GATA4, and GATA2. All have lower chromVAR scores in LOY. Related GATA motifs are correlated evidence and do not establish six independently identified regulators.</p>',
              figure('motifs/TF_extension/HOMER/known_motif_enrichment.png','Known-motif enrichment in gained and lost peaks using matched non-DA consensus backgrounds. One representative motif is shown per mapped group; displayed statistics remain motif-level.'),
              figure('motifs/TF_extension/concordance/HOMER_chromVAR.png','Known-motif enrichment versus per-sample accessibility change. Gained and lost peak subsets are shown separately. Agreement is method robustness on shared data, not independent validation.')]
        else:
            pieces += ['<p>HOMER enrichment is pending; prepared peak/background files are not enrichment results.</p>']
        if (tf_dir/'TOBIAS/summary.json').exists():
            footprints=pd.read_csv(tf_dir/'TOBIAS/TF_scores_with_chromVAR.csv')
            pieces += ['<h3>TOBIAS footprinting and method comparison</h3>',
              '<p>TOBIAS 0.17.5 applies Tn5 bias correction, footprint scoring, and BINDetect to condition pools with approximately equal read contributions from each replicate. Pool counts and reference lengths are verified. Signals use the common tested consensus, blacklist, and JASPAR2024 motif definitions; chrY and mitochondrial sequences are excluded.</p>',
              table(pd.read_csv(tf_dir/'TOBIAS/replicate_pool_accounting.csv')),
              table(pd.read_csv(tf_dir/'TOBIAS/ATACorrect_summary.csv')[['condition','input_reads','cut_site_FRiP','post_to_pre_bias_variance_forward','post_to_pre_bias_variance_reverse']]),
              '<p>These are TOBIAS internal QC metrics. Cut-site FRiP differs from the fragment-overlap FRiP reported earlier. A lower post/pre bias-variance ratio indicates reduced sequence bias in the tool\'s diagnostic; it is not experimental validation of TF occupancy.</p>',
              '<p>All six motif models shared by significant primary HOMER and chromVAR results also have negative TOBIAS changes, consistent with lower predicted occupancy in LOY. TFAP2C motifs also show negative TOBIAS changes. Very small pooled-site P-values can accompany tiny effects: Zfx has a change of -0.00909 and native P=4.83e-17, while its independent chromVAR FDR is 0.925 and linked ZFX RNA FDR is 0.970. Differential accessibility and chromVAR provide the biological-replicate tests; pooled-site significance alone is not a biological hit criterion.</p>',
              figure('motifs/TF_extension/TOBIAS/TOBIAS_differential_volcano.png','TOBIAS motif-score comparison, LOY versus CTRL. P-values and the added BH correction describe pooled motif-site patterns and are not biological-replicate differential tests.'),
              figure('motifs/TF_extension/concordance/chromVAR_TOBIAS.png','chromVAR accessibility changes versus TOBIAS predicted occupancy scores. Representatives are selected by chromVAR FDR within mapped groups before examining agreement.'),
              figure('motifs/TF_extension/TOBIAS/CTCF_aggregate.png','CTCF motif-site aggregate signal before and after Tn5 correction, using the same motif sites in both conditions.'),
              table(footprints.sort_values('pooled_site_BH_FDR')[['motif_id','motif_name','LOY_CTRL_change','LOY_CTRL_pvalue','pooled_site_BH_FDR','LOY_minus_CTRL_zscore','adj.P.Val']].head(20)),
              '<p>Footprints are occupancy predictions; they do not establish experimental binding or identify a unique paralog from shared sequence motifs. Pooling does not add biological replication. Full results and group-level TF RNA links are saved under motifs/TF_extension/.</p>']
        else:
            pieces += ['<p>TOBIAS footprinting is pending. Pooled footprint statistics will be interpreted as predicted occupancy and will not substitute for replicate-aware differential tests.</p>']
    pieces += ['<h2>RNAseq and lpWGS integration</h2>']
    rna_summary_path = ROOT/'integration/RNAseq/summary.json'
    integration_status = 'pending dataset selection and matching'
    if rna_summary_path.exists():
        rna_summary = json.loads(rna_summary_path.read_text())
        joint = pd.read_csv(ROOT/'integration/RNAseq/ATAC_RNA_joint_significant_genes.csv')
        integration_status = 'RNA group-effect integration complete; lpWGS and direct sample matching pending'
        pieces += [
          '<p>Selected RNA result: HT1080 combined-batch DESeq2 <code>~ batch + treatment</code>, LOY versus CTRL, following the user\'s instruction to use the result from Ben\'s thesis. The source notebook verifies this model and contrast; the thesis manuscript itself was not located.</p>',
          f'<p class="muted">Source: {escape(rna_summary["rna_source"])}<br>SHA256: {escape(rna_summary["rna_source_sha256"])}</p>',
          '<p>The original RNA model includes six CTRL and six LOY samples across AUG_2025 and MAY_2026, plus three parental samples contributing to model estimation. Existing RNA statistics were reused without refitting. This table uses EnsDb v86 annotation; the earlier voom v113 candidate was superseded. The source column named logCPM contains DESeq2 baseMean and is relabeled in the integration output.</p>',
          f'<p><b>{rna_summary["common_eligible_genes"]:,} genes</b> are evaluable in both assays. <b>{rna_summary["joint_significant_genes"]:,}</b> have RNA FDR &lt;0.05 and at least one promoter peak with ATAC FDR &lt;0.05: {rna_summary["concordant_gain"]} concordant gains, {rna_summary["concordant_loss"]} concordant losses, {rna_summary["discordant"]} discordant, and {rna_summary["mixed_promoter"]} mixed promoter patterns. Applying the RNA notebook\'s additional |log2FC| &gt;0.75 subset retains {rna_summary["joint_genes_RNA_abs_logFC_gt075"]} genes.</p>',
          '<p>Joining uses existing Entrez IDs and chromosome agreement. Ambiguous RNA mappings, missing statistics, and chrY/PAR genes are accounted for explicitly. The ATAC gene effect is a descriptive median across tested promoter peaks; minimum peak FDR is not a gene-level FDR. Concordance uses significant-promoter direction and RNA direction. All promoter links and separate exploratory nearest-TSS links are saved.</p>',
          table(pd.read_csv(ROOT/'integration/RNAseq/mapping_summary.csv')),
          figure('integration/RNAseq/ATAC_RNA_effects.png','Group-effect comparison on the common eligible gene universe. Colored genes have significant RNA and at least one significant ATAC promoter peak. Correlation and overlap counts are descriptive.'),
          table(joint[['rna_symbol','entrez_id','atac_median_promoter_logFC','atac_minimum_peak_FDR','rna_logFC','rna_FDR','joint_category']].head(25)),
          '<p>Full tables, the source manifest, mapping audit, and integration verification are in integration/RNAseq/. Direct ATAC/RNA sample matching is not established, so no sample-wise correlations are reported.</p>']
    else:
        pieces += ['<p>RNA dataset selection and sample matching remain pending.</p>']
    feature_dir = ROOT/'integration/feature_overlaps'
    if (feature_dir/'matched_FDR_overlap_grid.png').exists():
        pieces += ['<h2>Genomic feature overlaps by FDR</h2>',
          '<p>These diagrams compare distinct genes linked to significant ATAC peaks with significant RNA genes. Promoter peaks retain all promoter-overlap gene links, gene-body peaks retain all overlapping gene bodies, and intergenic peaks use exploratory nearest-TSS assignments. Each row uses its own common eligible gene universe; genes can occur in more than one feature row.</p>',
          figure('integration/feature_overlaps/matched_FDR_overlap_grid.png','Both ATAC and RNA use the displayed FDR cutoff: 0.05, 0.10, or 0.20. Any effect direction is included; no fold-change cutoff is imposed. Circle areas are schematic.'),
          figure('integration/feature_overlaps/RNA_fixed_005_overlap_grid.png','RNA stays at FDR <0.05 while ATAC is relaxed. RNA-only genes are restricted to the feature-specific eligible universe. ATAC significance remains peak-level; the overlap does not have a calibrated joint gene-level FDR.'),
          '<p>Individual feature plots, PNG and PDF exports, overlapping gene lists, mapping accounting, and the 18-panel count table are saved in integration/feature_overlaps/.</p>']
    scatter_dir = ROOT/'integration/feature_scatters'
    if (scatter_dir/'matched_FDR_scatter_grid.png').exists():
        pieces += ['<h2>Genomic feature ATAC–RNA scatter plots</h2>',
          '<p>Each point is one eligible gene from the corresponding overlap diagram. The x-axis is the median ATAC log2 fold change across all tested linked peaks in that feature; the y-axis is the existing batch-adjusted RNA log2 fold change. Colors show ATAC-only, RNA-only, both, or neither. All point coordinates and axis limits stay fixed across thresholds.</p>',
          figure('integration/feature_scatters/matched_FDR_scatter_grid.png','Promoter, gene-body, and intergenic gene effects with both assays at FDR <0.05, <0.10, and <0.20. Purple points are the overlapping gene sets from the Venn diagrams.'),
          figure('integration/feature_scatters/RNA_fixed_005_scatter_grid.png','The same feature-specific genes, with RNA FDR held at 0.05 while ATAC is relaxed. All 18 scatter-panel set counts exactly reproduce the overlap figures.'),
          '<p>ATAC significance refers to individual peaks; a median ATAC gene effect has no gene-level FDR and can differ in direction from an individual significant peak. Correlations use all eligible genes in each feature and are descriptive. Intergenic nearest-TSS links remain exploratory. Individual feature PNG/PDF plots, gene coordinates, memberships, and verification are in integration/feature_scatters/.</p>']
    review_dir = ROOT/'integration/promoter_FDR010_review'
    if (review_dir/'summary.json').exists():
        review = json.loads((review_dir/'summary.json').read_text())
        genes = pd.read_csv(review_dir/'all_94_genes.csv')
        categories = pd.read_csv(review_dir/'category_summary.csv')
        go = pd.read_csv(review_dir/'GO_BP_all_terms.csv')
        chr19_down = ((genes.status == 'Down in both') & (genes.chromosome == 'chr19')).sum()
        pieces += ['<h2 id="promoter-010-review">Promoter genes up and down at FDR 0.10/0.10</h2>',
          f'<p>Among {review["total"]} overlapping genes, <b>{review["up_both"]} increase in both assays</b>, <b>{review["down_both"]} decrease in both</b>, and <b>{review["opposing"]} have opposing directions</b>. By RNA direction alone, {review["RNA_up"]} are up and {review["RNA_down"]} are down. All effects compare LOY with CTRL using the selected batch + treatment RNA result.</p>',
          table(categories),
          '<p>The decreased set is more stable to tightening the threshold: 44 of 76 also pass 0.05/0.05, compared with 1 of 14 increased genes (LAYN). No minimum fold-change filter is imposed. The extra genes at 0.10 are exploratory candidates.</p>']
        for category, label in [('concordant_gain', 'Up in both'), ('concordant_loss', 'Down in both'), ('discordant', 'Opposing directions')]:
            selected = genes.loc[genes.category_at_cutoff == category]
            pieces += [f'<h3>{label}</h3>',
              figure(f'integration/promoter_FDR010_review/{category}_effects.png',
                'Paired ATAC promoter-median and RNA log2 fold changes, ordered by RNA FDR. The ATAC median is descriptive; significance is assessed at individual peaks.'),
              f'<details><summary>Gene statistics ({len(selected)} genes)</summary>',
              table(selected[['rna_symbol','chromosome','gene_description','atac_median_promoter_logFC','atac_minimum_peak_FDR','rna_logFC','rna_FDR','core_FDR005_both']]),
              '</details>']
        pieces += [f'<p>{chr19_down} of the 76 decreased genes are on chromosome 19. The decreased set contains 17 genes with ZNF-prefixed symbols. This positional pattern makes lpWGS copy-number context useful; it does not establish a copy-number explanation.</p>',
          f'<p>Exploratory GO biological-process enrichment returned <b>{review["GO_terms_FDR005"]} terms at FDR &lt;0.05</b>; the lowest BH-adjusted value across up and down tests together is {go.BH_FDR_across_up_down.min():.3f}. Guidance-related terms involving RAC3, GDNF, and DCC rank near the top of the increased set, but do not pass correction. These results do not establish pathway activation.</p>',
          '<details><summary>Top exploratory GO terms</summary>',
          table(go[['category','Description','Count','BH_FDR_across_up_down','query_symbols']].head(12)), '</details>',
          '<p>GO uses the common 11,898-gene assay-eligible universe (10,105 with GO annotation), term sizes 10-500, and BH correction across returned tests from both directions. Promoter-peak-number bias, annotation coverage, and related GO terms limit interpretation.</p>',
          '<p>Full gene lists, individual PNG/PDF plots, all linked promoter peaks, GO tables, and methods are in <code>integration/promoter_FDR010_review/</code>.</p>']
    wgs_status = '<p>lpWGS dataset selection and assay matching remain pending. No copy-number adjustment has been applied. Candidate hg38 WGS inputs are recorded in PROCESSING_PLAN.md; T2T WGS coordinates require conversion or recomputation before overlap.</p>'
    if (ROOT/'report/RNA_interaction_CN.html').exists():
        wgs_status = '<p>hg38 D0 WGS copy-number context and the formal RNA batch interaction are complete. <a href="RNA_interaction_CN.html">Open the interaction and copy-number report</a>. Direct culture matching remains unconfirmed; no sample-specific copy-number correction has been applied.</p>'
    pieces += [wgs_status,
      '<h2>Reproducibility and processing choices</h2>',
      '<p>Source BAM SHA256 hashes, BAM dictionaries, sample identities, reference checksums, software versions, and exact MACS3 commands are saved in config/, qc/, and logs/. Processing validates complete proper pairs, requires MAPQ ≥30 in both mates, excludes failed-QC/duplicate/nonprimary alignments, noncanonical contigs, and whole fragments overlapping the blacklist. Original alignments are not shifted or realigned.</p>',
      '<p>MACS3 3.0.4 uses BAMPE, q=0.01, explicit genome size 2,913,022,398, and the previously deduplicated fragments. Consensus candidates are 501 bp around summits, with summit support in at least two samples within either condition. Candidate ranking uses within-sample significance percentiles and deterministic greedy nonoverlap. Counts are integer fragment overlaps per feature; a fragment may overlap multiple features, while FRiP counts it once across the union. chrY and pseudoautosomal peaks remain available separately and are excluded from primary differential testing.</p>',
      '<p>Methods references: <a href="https://macs3-project.github.io/MACS/docs/callpeak.html">MACS3</a>, '
      '<a href="https://www.encodeproject.org/atac-seq/">ENCODE ATAC QC</a>, '
      '<a href="https://bioconductor.org/packages/release/bioc/html/csaw.html">csaw</a>, '
      '<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC9701614/">ATAC normalization benchmark</a>, '
      '<a href="https://mart.ensembl.org/info/genome/genebuild/human_PARS.html">GRCh38 PAR coordinates</a>.</p>',
      '</body></html>']
    if (ROOT/'report/TFAP2C_review.html').exists():
        pieces.insert(-1, '<section id="tfap2c-focused-review"><h2>Focused TFAP2C review</h2><p>All three TFAP2C motif models support selective AP-2-associated accessibility losses. Primary chromVAR FDR is 0.058–0.067; pooled predicted bound-site counts decrease by 2–6%. TFAP2C RNA is lower, while its annotated promoter peaks are not significantly less accessible. <a href="TFAP2C_review.html">Open the detailed cross-method review, replicate plots, footprint profiles and candidate genes.</a></p></section>')
    enlarged = 'integration/feature_scatters/promoter_ATAC010_RNA005_enlarged'
    if (ROOT/(enlarged+'.png')).exists():
        pieces.insert(-1, '<section id="promoter-atac010-rna005"><h2>Promoter overlap: ATAC FDR 0.10 / RNA FDR 0.05</h2><p>Enlarged selected panel: 11,898 eligible genes and 82 overlapping genes; no fold-change filter.</p>' + figure(enlarged+'.png', 'Original gene coordinates and full axis limits retained. ATAC membership uses any linked promoter peak; its median effect has no gene-level FDR.') + f'<p><a href="../{enlarged}.png">High-resolution PNG</a> · <a href="../{enlarged}.pdf">Vector PDF</a> · <a href="../{enlarged}_overlap_genes.csv">82 overlap genes</a></p></section>')
    publication = 'integration/feature_scatters/publication/promoter_ATAC010_RNA005_labeled'
    if (ROOT/(publication+'.png')).exists():
        pieces.insert(-1, '<section id="promoter-publication"><h2>Publication figure: promoter accessibility and RNA</h2><p>ATAC FDR &lt;0.10 and RNA FDR &lt;0.05: 82 overlapping genes (9 concordant increases, 70 concordant decreases, 3 opposing directions). Twenty-four labels show 12 concordant overlap genes in dark text plus six ATAC-only and six RNA-only outliers in lighter blue/orange.</p>' + figure(publication+'.png', 'All eligible gene coordinates retained. Labels are illustrative selections; promoter median effects have no gene-level FDR.') + f'<p><a href="../{publication}.pdf">Vector PDF</a> · <a href="../{publication}.svg">SVG</a> · <a href="../{publication}.png">600 dpi PNG</a> · <a href="../integration/feature_scatters/publication/figure_caption.txt">Figure legend</a></p></section>')
    if (ROOT/'report/down70_enrichment.html').exists():
        pieces.insert(-1, '<section id="down70-enrichment"><h2>Enrichment of the 70 concordant down/down genes</h2><p>No GO term passes FDR 0.05. Chromosome 19 is enriched (19 genes; 3.92-fold; FDR 3.72e-6), and the ZNF/ZFP-named component remains enriched after chromosome conditioning. Local 1 Mb clustering is significant globally (FDR 0.00330), with shared promoter peaks limiting independence. <a href="down70_enrichment.html">Full functional and chromosome analysis</a>.</p></section>')
    if (ROOT/'report/batch_RNA_overlaps.html').exists():
        pieces.insert(-1, '<section id="batch-rna-overlaps"><h2>ATAC overlap with each RNA-seq batch separately</h2><p>Separate August 2025 and May 2026 DESeq2 fits, ATAC FDR &lt;0.10 and RNA FDR &lt;0.05. <a href="batch_RNA_overlaps.html">Overlap plots, scatter plots, shared genes and support for the original 70 down/down genes</a>.</p></section>')
    if (ROOT/'report/gene_atlas.html').exists():
        pieces.insert(-1, '<section id="gene-locus-atlas"><h2>ATAC gene-locus diagrams</h2><p>Gene overviews and peak close-ups for strong DAR-linked genes and ATAC/RNA overlap genes at FDR 0.10/0.05. <a href="gene_atlas.html">Search the atlas, inspect all replicates, and download diagrams and peak measurements</a>.</p></section>')
    if (ROOT/'report/interactive_ATAC_RNA.html').exists():
        pieces.insert(-1, '<section id="interactive-atac-rna"><h2>Interactive ATAC versus RNA explorer</h2><p><a href="interactive_ATAC_RNA.html">Open the interactive scatter, sample expression and ATAC tracks</a>. Click any eligible gene, switch RNA datasets or genomic features, and inspect individual peaks. For on-demand coverage, start the local viewer with <code>Start_interactive_report.cmd</code> in the processing folder.</p></section>')
    if (ROOT/'report/pathway_GSEA.html').exists():
        pieces.insert(-1, '<section id="pathway-gsea"><h2>ATAC and RNA pathway enrichment</h2><p>Hallmark and Reactome GSEA on the common promoter-linked gene universe. <a href="pathway_GSEA.html">Pathway heatmaps, ATAC–RNA enrichment comparisons, shared leading-edge genes and sensitivity analyses</a>.</p></section>')
    if (ROOT/'report/pathway_gene_drivers.html').exists():
        pieces.insert(-1, '<section id="pathway-gene-drivers"><h2>Inspect the genes driving enrichment</h2><p><a href="pathway_gene_drivers.html">E2F, MYC and interferon: gene effects, individual RNA and ATAC samples, and interactive promoter peak inspection</a>.</p></section>')
    (ROOT/'report/ATACseq_report.html').write_text('\n'.join(pieces),encoding='utf-8')
    overview = {'consensus_regions':len(consensus),'tested_regions':int(primary_summary.tested),
      'primary_gained_FDR005':int(primary_summary.gained),'primary_lost_FDR005':int(primary_summary.lost),
      'integration_status':integration_status}
    (ROOT/'report/summary.json').write_text(json.dumps(overview,indent=2)+'\n')
    print(json.dumps(overview,indent=2))
