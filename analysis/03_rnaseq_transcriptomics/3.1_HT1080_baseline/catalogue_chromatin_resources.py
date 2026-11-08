"""Build curated resource tables from the completed read-only discovery audit."""
from pathlib import Path
import hashlib
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'provenance/chipseq_resources'


def main():
    metadata = pd.read_csv(OUT/'geo_local_sample_metadata.csv').fillna('')
    inventory = pd.read_csv(OUT/'upstream_file_inventory.csv')
    rows = []
    selected = metadata[metadata.metadata_source.str.contains('GSE261834') &
                        metadata.source_name_ch1.eq('HT1080')]
    for _, sample in selected.iterrows():
        files = inventory[inventory.relative_path.str.contains(sample.accession, regex=False)]
        rows.append({
            'accession': sample.accession, 'title': sample.title,
            'biosample': sample.source_name_ch1,
            'characteristics': sample.characteristics_ch1,
            'assembly': 'hg38', 'assay': sample.library_strategy,
            'local_file_count': len(files),
            'local_files': ';'.join(files.absolute_path),
            'metadata_source': sample.metadata_source,
            'public_url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sample.accession}',
            'caution': ('Title/condition discrepancy: condition says H3K4me4; title identifies mark. '
                        'Confirm antibody/source supplement before use.' if sample.accession in
                        ['GSM8152733', 'GSM8152737'] else
                        'Processed files may pool replicates; file accession alone does not resolve pooling.')
        })
    pd.DataFrame(rows).to_csv(OUT/'GSE261834_HT1080_sample_manifest.csv', index=False)

    mark_rows = []
    for mark, first, reps, form in [
        ('H2AK119ub1',8152720,4,'broad'), ('H3K27me3',8152724,2,'broad'),
        ('H3K9me3',8152726,2,'broad'), ('H4K20me3',8152728,2,'broad'),
        ('H4K20me1',8152730,2,'broad'), ('H3K4me1',8152732,2,'narrow'),
        ('H3K27ac',8152734,2,'narrow'), ('H3K4me3',8152736,2,'narrow')]:
        acc = f'GSM{first}'
        subset = inventory[inventory.relative_path.str.contains(f'GSE261834_RAW/{acc}_', regex=False)]
        mark_rows.append({'mark': mark, 'source_accessions': ';'.join(f'GSM{i}' for i in range(first,first+reps)),
                          'metadata_replicates': reps, 'processed_accession': acc,
                          'bigwig_local': any(subset.relative_path.str.endswith('.bw')),
                          'peak_format': form,
                          'peak_file_local': any(subset.relative_path.str.contains('Peak')),
                          'replicate_pooling': 'Unresolved; only first accession has processed tracks locally',
                          'files': ';'.join(subset.absolute_path)})
    pd.DataFrame(mark_rows).to_csv(OUT/'GSE261834_measured_mark_panel.csv', index=False)

    encode = []
    for path in sorted((OUT/'public_metadata').glob('ENCFF*.json')):
        d = json.loads(path.read_text(encoding='utf-8'))
        encode.append({'accession': d['accession'], 'biosample': d.get('biosample_ontology',{}).get('term_name'),
                       'assembly': d.get('assembly'), 'annotation_type': d.get('annotation_type'),
                       'output_type': d.get('output_type'),
                       'target': ';'.join(t.get('label','') for t in d.get('targets',[])),
                       'submitted_file_name': d.get('submitted_file_name'),
                       'dataset': d.get('dataset'), 'local_metadata': str(path),
                       'public_url': f'https://www.encodeproject.org/files/{d["accession"]}/'})
    pd.DataFrame(encode).to_csv(OUT/'ENCODE_provenance_audit.csv', index=False)

    studies = [
      ('GSE261834','HT1080; subset excludes K562','Measured reference','Priority 1',
       'NEED-seq histone panel, lamins, TFs; EM-seq and Hi-C',
       '80 files under source dataset, including archive; 8 standard histone marks have peaks and bigWigs',
       'hg38','Replicate pooling and assay controls need resolution; EM-seq calls absent locally',
       'Build primary promoter/region annotations independently of RNA results',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261834'),
      ('GSE135580','HT1080 subset; also HCT116 and DKO','Measured reference','Priority 2',
       'H3K27me3/H3K9me2/H3K9me3/HP1alpha/CTCF ChIP; Protect-seq; Repli-seq',
       'GEO metadata and rendered report local; signal tracks referenced on HPC not found in inspected local data folder',
       'hg38','Some processed domain tracks are coarse 25-kb bins; inspect signal resolution',
       'Cross-study support for Polycomb versus constitutive heterochromatin',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135580'),
      ('GSE153869','HT1080 WT and transfected derivatives','Measured reference and perturbation','Priority 2',
       'H3K27me3/H3K27ac/CTCF/DDIT3 ChIP',
       'Public BED/raw data; not in inspected ChIPseq data folder',
       'hg38','Subset WT controls; do not pool DDIT3, FUS-DDIT3, or EWS-FLI1 conditions',
       'Independent WT H3K27me3 and H3K27ac validation',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153869'),
      ('GSE86504','HT1080','Measured reference and perturbation','Priority 3',
       'H3K4me3, MYC, MEN1, Pol II ChIP; MYC/MEN1 knockdown RNA',
       'Public BED/raw data; not in inspected ChIPseq data folder',
       'hg19','Resolve generic shRNA metadata for baseline ChIP; harmonize assembly',
       'Validate active promoters and explore MYC/MEN1 regulation',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86504'),
      ('GSE162690','HT1080 GSM4957247/48 among 10 cell lines','Measured accessibility reference','Priority 3',
       'Bulk ATAC duplicates and mixed-cell scATAC',
       'Metadata and old report local; HPC tracks not present in inspected local folder',
       'hg19 in GEO','Old report uses EnsDb v86 hg38 TSS; original HPC build conversion unverified',
       'Independent basal accessible regions after build verification',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162690'),
      ('GSE181247','HT1080','Measured perturbation','Priority 3',
       'ATAC across 3 collagen concentration/confinement conditions, 4 replicates each',
       'Metadata and rendered report local; HPC tracks not present in inspected local folder',
       'Verify before overlap','Confinement is an experimental perturbation, not a LOY comparison',
       'Compare ECM/migration accessibility programs with current RNA directions',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181247'),
      ('ENCSR295VRM / ENCFF461DZD','HT1080 BSS00709','Model-derived chromatin states','Secondary sensitivity',
       'EpiMap ChromHMM 18 states',
       'Verified public metadata cached; segmentation itself not downloaded',
       'GRCh38 available; also hg19','EpiMap combines observed and imputed inputs; per-mark observed status unresolved',
       'Compare promoter annotations with directly measured NEED-seq marks',
       'https://www.encodeproject.org/annotations/ENCSR295VRM/'),
      ('ENCODE Avocado tracks cited in gkaf720','HT1080','Imputed predictions','Secondary sensitivity only',
       'H3K4me3, H3K27ac, H3K27me3, H3K9ac, ATAC and POLR2A predictions',
       'Six file metadata records verified and cached',
       'GRCh38','ENCFF408PEC is H3K9ac, although paper calls it H3K9me3; all six are imputed',
       'Context only, never independent experimental confirmation',
       'https://www.encodeproject.org/files/ENCFF633FBI/'),
      ('10.64898/2026.04.08.717275','HT1080','Measured perturbation described in preprint','High-interest acquisition lead',
       'HDAC/LSD1 inhibition; NicE-seq; NEED-seq H3/active marks/RCOR1/RUNX/JunB; RNA',
       'Full-text cached; no GEO accession supplied in inspected full text',
       'Not frozen','Preprint; data location unresolved; Results/Methods differ on LSD1 exposure duration',
       'External HT1080 HDAC/CoREST response benchmark if deposited data are located',
       'https://pmc.ncbi.nlm.nih.gov/articles/PMC13081975/'),
      ('GSE82260','DU145','External mechanistic context','Optional',
       'ZMYND8/JARID1D(KDM5D), histone marks and knockdown RNA',
       'Public; not an HT1080 reference',
       'Verify if acquired','Different cancer lineage; no transfer of binding identity without evidence',
       'Contrast mechanism-derived targets with directly measured HT1080 promoter states',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82260'),
      ('GSE147672','Human brain','Cross-tissue comparison','Context only',
       'Brain accessibility/3D regulatory atlas',
       'Metadata and rendered brain comparison local',
       'Verify original track build','Not HT1080',
       'Optional neural-program context; exclude from HT1080-specific sets',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147672'),
      ('GSE324952','Mouse 10T1/2','Cross-species perturbation','Context only',
       'Ranbp2 wild-type versus knockout ATAC',
       'Rendered mouse-ortholog comparison local; no source tracks in inspected data folder',
       'Mouse; verify exact assembly','Not human or HT1080',
       'Optional ortholog-based context; exclude from HT1080 reference',
       'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324952'),
      ('Active Motif Q54521 / SEPT_08_2026_REPROCESS','Own HT1080 CTRL and LOY','Measured LOY accessibility contrast','Group-effect integration; design revision needed',
       '6 ATAC libraries; fresh peaks, edgeR, csaw, HOMER, chromVAR, TOBIAS',
       'Original inputs and newer completed reanalysis in projects/ATACseq; multiple pipelines use same libraries',
       'hg38','User confirmed successive passages on 2026-09-14; direct RNA partners unknown. Prior independent-replicate ATAC inference and older RNA integration require revision',
       'Join group-level accessibility effects to current RNA and reference promoter features',
       'A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS/PROCESSING_PLAN.md'),
    ]
    columns = ['resource','biosample','evidence_type','priority','assays','availability',
               'assembly','limitations','proposed_use','source']
    pd.DataFrame(studies, columns=columns).to_csv(OUT/'chromatin_resource_catalogue.csv', index=False)
    # Hash the small source documents that support local provenance conclusions.
    source_files = [Path('A:/multiome/scripts/ChIPseq/config_GSE261834.R'),
                    Path('A:/multiome/scripts/ChIPseq/01_bivalent_loy_overlap.Rmd'),
                    Path('A:/multiome/results/ChIPseq/tables/bivalent_genes_HT1080.csv'),
                    Path('A:/multiome/results/ChIPseq/tables/de_HT1080_LOY_chromatin_annotated.csv'),
                    Path('A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS/PROCESSING_PLAN.md'),
                    Path('A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS/RUN_STATUS.md')]
    source_files += [Path(p) for p in metadata.metadata_source.unique()]
    hashes = [{'path': str(p), 'bytes': p.stat().st_size,
               'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in source_files]
    pd.DataFrame(hashes).to_csv(OUT/'audit_source_checksums.csv', index=False)
    # Remove the old overly broad regex match GSEA from the initial inventory.
    accs = pd.read_csv(OUT/'legacy_report_accessions.csv')
    accs[accs.accession.str.match(r'(?:(?:GSE|GSM|SRR|SRX)\d+|(?:ENCSR|ENCFF)[A-Z0-9]+)$')].to_csv(
        OUT/'legacy_report_accessions.csv', index=False)
    print(json.dumps({'HT1080_NEED_series_samples':len(rows),'standard_marks':len(mark_rows),
                      'curated_resources':len(studies),'encode_records':len(encode)}))


if __name__ == '__main__':
    main()
