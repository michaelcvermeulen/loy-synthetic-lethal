$ErrorActionPreference = 'Stop'

$analysis = Split-Path -Parent $PSScriptRoot
$thesis = 'A:\Dropbox\LOY_paralogs\projects\SEPT_11_2026_THESIS'
$utf8 = New-Object System.Text.UTF8Encoding($false)

function Read-Utf8([string]$path) {
    return [System.IO.File]::ReadAllText($path, [System.Text.Encoding]::UTF8)
}

function Write-Utf8([string]$path, [string]$text) {
    [System.IO.File]::WriteAllText($path, $text, $utf8)
}

function Replace-Exact([string]$path, [string]$old, [string]$new) {
    $text = Read-Utf8 $path
    $oldWindows = $old.Replace("`r`n", "`n").Replace("`n", "`r`n")
    $newWindows = $new.Replace("`r`n", "`n").Replace("`n", "`r`n")
    if ($text.Contains($new) -or $text.Contains($newWindows)) {
        return
    } elseif ($text.Contains($old)) {
        Write-Utf8 $path ($text.Replace($old, $newWindows))
    } elseif ($text.Contains($oldWindows)) {
        Write-Utf8 $path ($text.Replace($oldWindows, $newWindows))
    } else {
        throw "Expected text not found in $path"
    }
}

function Replace-ParagraphStarting([string]$path, [string]$prefix, [string]$new) {
    $text = Read-Utf8 $path
    if ($text.Contains($new)) { return }
    $start = $text.IndexOf($prefix)
    if ($start -lt 0) { throw "Expected paragraph prefix not found in $path" }
    $finish = $text.IndexOf("`r`n`r`n", $start)
    if ($finish -lt 0) { $finish = $text.Length }
    $newWindows = $new.Replace("`r`n", "`n").Replace("`n", "`r`n")
    Write-Utf8 $path ($text.Substring(0, $start) + $newWindows + $text.Substring($finish))
}

function Append-Once([string]$path, [string]$marker, [string]$addition) {
    $text = Read-Utf8 $path
    if ($text.Contains($marker)) { return }
    $additionWindows = $addition.Replace("`r`n", "`n").Replace("`n", "`r`n")
    Write-Utf8 $path ($text.TrimEnd() + "`r`n`r`n" + $additionWindows.Trim() + "`r`n")
}

$resultsPath = Join-Path $thesis 'Chapter3_results.md'
$methodsPath = Join-Path $thesis 'Chapter3_methods.md'
$captionsPath = Join-Path $thesis 'Chapter3_figure_captions.md'
$sourceNotesPath = Join-Path $thesis 'Chapter3_source_notes.md'
$resultsPlanPath = Join-Path $thesis 'Chapter3_results_plan.md'
$agentsPath = Join-Path $thesis 'AGENTS.md'

$oldResults = @'
### 3.4.8 Transcriptional responses to UTY and KDM5D re-expression

*[Establish construct identity and expression, vector controls and the Y-retaining reference; then describe reversal, persistence, overshoot and new changes after individual or combined re-expression. Separate molecular reversal from functional dependency rescue.]*

*The [molecular Results plan](research/model_development/genomic_molecular_characterisation_plan.md) specifies main and supplementary displays. Bracketed passages are drafting placeholders, not reported findings.*
'@
$newResults = (Read-Utf8 (Join-Path $analysis 'thesis_materials\results_3_4_8.md')).Trim()
if (-not (Read-Utf8 $resultsPath).Contains('We next asked whether re-expression of two Y-linked chromatin regulators altered the transcriptional state')) {
    Replace-Exact $resultsPath $oldResults.Trim() $newResults
}

$methods = Read-Utf8 (Join-Path $analysis 'thesis_materials\methods_3_3_15.md')
Append-Once $methodsPath '### 3.3.15 HT1080 UTY and KDM5D re-expression RNA analysis' $methods

$captions = Read-Utf8 (Join-Path $analysis 'thesis_materials\figure_captions.md')
Append-Once $captionsPath '## main - HT1080 add-back design, construct expression and global RNA structure' $captions

$oldPlan = '**Source and readiness:** S28 is the designated final source, supported by S03/S04/S21 and source-matched S23 synthesis.'
$newPlan = '**Source and readiness:** `RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026` is the designated final numerical source, supported by S03/S04/S21 and source-matched S23 synthesis. The frozen thesis analysis supersedes the earlier S28 `reannotate_v105` designation; S28 and the September consolidated pipeline remain provenance and sensitivity sources.'
Replace-Exact $resultsPlanPath $oldPlan $newPlan

$newSource = '**Status: `EXISTS / RECONCILED / WRITTEN`.** Use `HT1080_ADDBACK/THESIS_ADDBACK_2026` for the Chapter 3 numerical add-back results. It uses the corrected sample identities, GENCODE v49 annotation, all 15 pairwise contrasts and the August-only frozen baseline signature. The prior `reannotate_v105` notebooks, September consolidated analysis and review reports remain provenance and sensitivity sources rather than final numerical sources. See [the add-back source notes](research/molecular_characterisation/ht1080_addback_results_notes.md).'
Replace-ParagraphStarting $sourceNotesPath '**Status: `EXISTS / RECONCILE / WRITE`.** Use the user-designated final `HT1080_ADDBACK/reannotate_v105`' $newSource

$oldFinish = '**Output:** Figure 3.11, Table 3.6, Figure S3.15 and Table S3.15. **Sources:** S03/S04 for biological units and provenance, S21 for historical label/regional-Y context, S28 for final add-back sources and S23 for qualified molecular synthesis. **Finish:** freeze construct mapping, actual protein-validation status, primary contrast/signature definitions and exact effect tables. No re-expression experiment or RNA model is rerun for this planning update.'
$newFinish = '**Output:** Two main add-back figures, one supplementary sensitivity figure and complete source tables. **Sources:** S03/S04 for biological units and provenance, S21 for historical label/regional-Y context, the frozen thesis analysis for numerical results and S23 for qualified molecular synthesis. **Remaining finish:** replace the reserved Western panel after checking original images, lane identities, loading controls and exposure information.'
Replace-Exact $sourceNotesPath $oldFinish $newFinish

$agentNote = @'
# HT1080 add-back thesis source

Author-confirmed 21 September 2026: the empty-vector, UTY, KDM5D and combined populations were produced by separate transductions from the same starting Y-depleted pool. Sample indices 1 to 3 are matched collection occasions across the six groups, not independent transductions. Western blotting was completed for UTY and KDM5D in each sample and showed the expected abundance; source images and lane metadata still require checking before final figure insertion. Use `RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026` as the final numerical source for Results 3.4.8. It supersedes the earlier `reannotate_v105` designation while preserving older analyses as provenance and sensitivity sources.
'@
Append-Once $agentsPath '# HT1080 add-back thesis source' $agentNote

# Repair an encoding artefact left by the first Windows PowerShell integration
# pass. The replacement is safe to rerun and preserves the rest of AGENTS.md.
$badDash = '1' + [char]0x00e2 + [char]0x20ac + [char]0x201c + '3'
Replace-Exact $agentsPath ("Sample indices " + $badDash + " are matched collection occasions") 'Sample indices 1 to 3 are matched collection occasions'

# Reconcile the older planning/source register with the implemented analysis.
# S28 remains useful provenance, while S31 records the frozen numerical source.
Replace-Exact $sourceNotesPath '### 3.3.18 HT1080 Y-gene add-back, construct verification and transcriptional reversal' '### 3.3.15 HT1080 Y-gene add-back, construct verification and transcriptional reversal'

$oldResultSource = @'
### 3.4.10 Re-expression of selected Y-linked genes tests the reversibility of the HT1080 transcriptional phenotype

**Status:** `EXISTS / RECONCILE / WRITE`. **Claim:** the add-back experiment assesses the extent and specificity of transcriptional reversal following UTY/KDM5D re-expression, with construct, vector and population context attached to that interpretation. The subsection should report the retained results without presupposing that either construct fully restores the control state.

**Exact source route:** S28 final main, construct-validation, plasmid-only and confounder notebooks and their output directories; S04/S21 label and regional-Y provenance; S23 qualified molecular synthesis where it uses the same final inputs. **Finish:** verify specimen/construct mapping and actual validation modalities, freeze the basal signature and primary add-back contrasts, extract exact effects and define reversal without equating it to functional rescue.
'@
$newResultSource = @'
### Add-back subsection implemented as Results 3.4.8

**Status:** `EXISTS / RECONCILED / WRITTEN`. **Claim:** UTY-containing established populations showed partial molecular reversal of a frozen LOY-associated expression signature, whereas KDM5D alone did not show a detectable corrective shift under the retained analyses. The inference is limited to differences among the sampled established populations.

**Exact source route:** S31 frozen thesis analysis for numerical results and figure source data; S03/S04 for biological units and provenance; S21/S28 for historical label, regional-Y and analysis context; S23 for qualified molecular synthesis. **Remaining finish:** replace the reserved Western-blot panel after checking the original images, lane identities, loading controls and exposure information.
'@
Replace-Exact $sourceNotesPath $oldResultSource.Trim() $newResultSource.Trim()

$oldFigureRow = '| **3.11** HT1080 Y-gene add-back | 3.4.10 | Construct verification, vector context, individual/dual effects and extent of transcriptional reversal | S28 with S03/S04/S21; freeze labels, contrasts, signatures and validation modalities. |'
$newFigureRow = '| **3.11** HT1080 Y-gene add-back | 3.4.8 | Construct verification, vector context, individual/dual effects and extent of transcriptional reversal | S31 with S03/S04/S21/S28; numerical analysis, labels and signature are frozen; Western source panel remains pending. |'
Replace-Exact $sourceNotesPath $oldFigureRow $newFigureRow

$oldS21 = '| **S21** | `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/`; `P/RNAseq/analysis/HT1080_ADDBACK/SEPT_2026_REANALYSIS/RUN_STATUS.md` | Source mapping for May basal libraries and regional-Y/label audit context. The final add-back branch is reannotate_v105, detailed in S28. Historical consolidated counts must not replace it. Molecular add-back now belongs in Chapter 3. |'
$newS21 = '| **S21** | `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/`; `P/RNAseq/analysis/HT1080_ADDBACK/SEPT_2026_REANALYSIS/RUN_STATUS.md` | Historical May-library mapping, regional-Y and label-audit context. The numerical Chapter 3 add-back source is S31; these branches remain provenance and sensitivity sources. |'
Replace-Exact $sourceNotesPath $oldS21 $newS21

$oldS23 = '| **S23** | `L/_thesis_synthesis/thesis_evidence_synthesis/SYNTHESIS.md`; `L/_thesis_synthesis/thesis_evidence_synthesis/ANALYSIS_PROTOCOL.md`; `L/_thesis_synthesis/thesis_evidence_synthesis/PROTOCOL_NOTES.md`; `L/_thesis_synthesis/thesis_evidence_synthesis/INPUTS.tsv`; `L/_thesis_synthesis/thesis_evidence_synthesis/results/` | Audited public A549 comparison (WP2) and molecular add-back synthesis (relevant WP1 components). Match selected programmes and input hashes to S28 before use. Keep screen-linked mechanism in Chapter 4; do not treat baseline ATAC overlap as add-back ATAC. |'
$newS23 = '| **S23** | `L/_thesis_synthesis/thesis_evidence_synthesis/SYNTHESIS.md`; `L/_thesis_synthesis/thesis_evidence_synthesis/ANALYSIS_PROTOCOL.md`; `L/_thesis_synthesis/thesis_evidence_synthesis/PROTOCOL_NOTES.md`; `L/_thesis_synthesis/thesis_evidence_synthesis/INPUTS.tsv`; `L/_thesis_synthesis/thesis_evidence_synthesis/results/` | Audited public A549 comparison (WP2) and historical molecular add-back synthesis (relevant WP1 components). Use S31 for final add-back numbers. Keep screen-linked mechanism in Chapter 4; do not treat baseline ATAC overlap as add-back ATAC. |'
Replace-Exact $sourceNotesPath $oldS23 $newS23

$oldS28 = '| **S28** | `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK_construct_validation.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK_plasmid_only.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK_confounders.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/DE/`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/DE_construct_validation/`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/figures_addback/`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/figures_construct_validation/` | Located final add-back source notebooks and outputs. Main and four-group subsets reuse observations and differ in model scope. Their historical claims about complete Y loss, pure vector effects or UTY catalytic function require qualification; use exact results with corrected provenance. |'
$newS28 = '| **S28** | `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK_construct_validation.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK_plasmid_only.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/analysis_HT1080_ADDBACK_confounders.Rmd`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/DE/`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/DE_construct_validation/`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/figures_addback/`; `P/RNAseq/analysis/HT1080_ADDBACK/reannotate_v105/figures_construct_validation/` | Legacy add-back notebooks and outputs retained for provenance and sensitivity. Main and four-group subsets reuse observations and differ in model scope. Historical claims about complete Y loss, pure vector effects or UTY catalytic function are not carried into the thesis interpretation. |'
Replace-Exact $sourceNotesPath $oldS28 $newS28

$s31 = '| **S31** | `P/RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/`; `P/SEPT_11_2026_THESIS/figures/chapter3/source_data/ht1080_addback_2026/` | Frozen final numerical source for Results 3.4.8: corrected labels, GENCODE v49 counts, collection-occasion-blocked contrasts, August-only frozen signature, sensitivity analyses, input hashes and figure source tables. |'
$s30Marker = '| **S30** | `L/code/find_unique_seq/`; linked `L/data/Y_seq_chm13.RDS`, `unique_Y_sequences_chm13.RDS`, `candidate_guides.csv`, `DeepSpCas9_Scores*.txt`, `top_hits.csv` | Original design directory supplied 16 September 2026. PDF/Rmd/code review and saved-artifact audit support expanded 3.4.5 and correct the 7,470 threshold/unit attribution. Exact files, roles and hashes: [source audit](research/sgRNA_design/sgRNA_design_source_audit_2026-09-16.md) and its manifest. No original genome scan or scoring pipeline was rerun. |'
Replace-Exact $sourceNotesPath $s30Marker ($s30Marker + "`r`n" + $s31)

$oldAnchor = '| HT1080 add-back sample layout | Source describes 18 libraries across six groups, including four vector-bearing groups | S28. These are library counts; trace population/construct units and the reuse of basal May samples separately. |'
$newAnchor = '| HT1080 add-back sample layout | 18 libraries across six groups; one established population per group sampled on three matched collection occasions | S31 with S03/S04. The occasions are not independent transductions; basal May samples are reused as contextual references. |'
Replace-Exact $sourceNotesPath $oldAnchor $newAnchor

$oldPriority = '| **1** | Molecular add-back source and definition | User-designated final reannotate_v105 sources exist | Freeze primary main/subset model, basal signature, direct construct contrasts and gene-level reversal rules with exact denominators; keep older rescue lists separately labelled. |'
$newPriority = '| **1** | Molecular add-back source and definition | S31 analysis, signature and contrasts are frozen and written into Results 3.4.8 | Preserve S31 as the final numerical source and keep older rescue lists separately labelled as provenance or sensitivity analyses. |'
Replace-Exact $sourceNotesPath $oldPriority $newPriority

$oldValidation = '| **1** | Construct labels, validation and controls | Dedicated expression checks and historical relabel provenance exist | Map original/corrected IDs, vector/selection exposures and pool relationships; report the actual RNA/protein evidence and the absent Y-intact vector comparison. |'
$newValidation = '| **1** | Construct labels, validation and controls | Corrected IDs, pool relationships and RNA/protein verification are documented | Insert the original Western images only after checking lane identities, loading controls and exposure information; retain the absence of a Y-intact vector control as a design limit. |'
Replace-Exact $sourceNotesPath $oldValidation $newValidation

$researchDir = Join-Path $thesis 'research\molecular_characterisation'
New-Item -ItemType Directory -Force -Path $researchDir | Out-Null
Copy-Item -LiteralPath (Join-Path $analysis 'thesis_materials\ht1080_addback_results_notes.md') -Destination (Join-Path $researchDir 'ht1080_addback_results_notes.md') -Force

$figDir = Join-Path $thesis 'figures\chapter3'
New-Item -ItemType Directory -Force -Path $figDir | Out-Null
Get-ChildItem -LiteralPath (Join-Path $analysis 'figures') -File | ForEach-Object {
    Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $figDir $_.Name) -Force
}

$sourceDir = Join-Path $figDir 'source_data\ht1080_addback_2026'
New-Item -ItemType Directory -Force -Path $sourceDir | Out-Null
$sourceFiles = @(
    'outputs\sample_manifest.csv',
    'outputs\filter_summary.csv',
    'outputs\de_summary.csv',
    'outputs\factorial_summary.csv',
    'outputs\construct_rna_abundance.csv',
    'outputs\pca_scores.csv',
    'outputs\pca_variance.csv',
    'outputs\claim_evidence.csv',
    'outputs\signature\frozen_signature_mapping.csv',
    'outputs\signature\frozen_signature_scores.csv',
    'outputs\signature\score_contrasts.csv',
    'outputs\signature\fry_directional.csv',
    'outputs\signature\distance_to_control_by_gene.csv',
    'outputs\signature\distance_to_control_summary.csv',
    'outputs\sensitivity\unblocked_summary.csv',
    'outputs\sensitivity\leave_one_occasion_out.csv',
    'provenance\input_manifest.csv',
    'provenance\verification.txt',
    'provenance\sessionInfo.txt',
    'README.md',
    'code\run_analysis.R'
)
foreach ($relative in $sourceFiles) {
    $from = Join-Path $analysis $relative
    $name = $relative.Replace('\', '__')
    Copy-Item -LiteralPath $from -Destination (Join-Path $sourceDir $name) -Force
}

Write-Output 'Thesis add-back section, methods, captions, notes, figures and source data integrated.'
