param([Parameter(Mandatory=$true)][ValidateSet('HT1080','DU145')][string]$Cell)
$ErrorActionPreference='Stop'
$analysisFolder='SEPT_11_2026_THESIS/research/model_development/ht1080_du145_ypar_20261027'
$samples = if ($Cell -eq 'HT1080') { @('HT1080-sgCTRL-D0','HT1080-LOY-D0') } else { @('LNCaP-sgCTRL','LNCaP-LOY-2X') }
$deadline=(Get-Date).AddHours(4)
while ($true) {
    $ready=@($samples | Where-Object { Test-Path -LiteralPath "WGS/reanalysis/Y_PAR_masked_20261027/qc/$_.REALIGNMENT_COMPLETE.txt" })
    if ($ready.Count -eq 2) { break }
    if ((Get-Date) -gt $deadline) { throw "Timed out waiting for $Cell validated alignments" }
    Start-Sleep -Seconds 5
}
Write-Output "$Cell alignments validated; starting matched-pair analysis."
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla "$analysisFolder/fit_pairs.R" $Cell masked *> "$analysisFolder/$Cell/masked_fit.log"
if ($LASTEXITCODE -ne 0) { throw "$Cell fit failed; see masked_fit.log" }
Write-Output "$Cell masked pair analysis complete."
