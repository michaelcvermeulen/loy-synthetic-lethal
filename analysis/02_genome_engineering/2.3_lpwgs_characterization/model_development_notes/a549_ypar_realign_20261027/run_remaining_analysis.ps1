$ErrorActionPreference = 'Stop'
$analysisFolder = 'SEPT_11_2026_THESIS/research/model_development/a549_ypar_realign_20261027'
$completeMarker = 'WGS/reanalysis/Y_PAR_masked_20261027/qc/A549-LOY.REALIGNMENT_COMPLETE.txt'
$deadline = (Get-Date).AddHours(2)
while (-not (Test-Path -LiteralPath $completeMarker)) {
    if ((Get-Date) -gt $deadline) { throw 'Timed out waiting for validated LOY alignment.' }
    Start-Sleep -Seconds 5
}
Write-Output 'Validated LOY alignment ready; counting reads.'
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla "$analysisFolder/count_new_bams.R" 'A549-LOY' *> "$analysisFolder/count_LOY.log"
if ($LASTEXITCODE -ne 0) { throw 'LOY read counting failed; see count_LOY.log.' }
Write-Output 'Counts complete; fitting the matched A549 pair.'
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla "$analysisFolder/fit_new_pair.R" *> "$analysisFolder/fit.log"
if ($LASTEXITCODE -ne 0) { throw 'Matched-pair analysis failed; see fit.log.' }
Write-Output 'Matched-pair QDNAseq analysis complete.'
