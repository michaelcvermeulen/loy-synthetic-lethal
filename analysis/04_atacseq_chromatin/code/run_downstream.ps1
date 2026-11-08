param([switch]$FinalizeOnly)
$ErrorActionPreference = 'Stop'
$atacRunRoot = Split-Path -Parent $PSScriptRoot
$atacRscript = 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe'
$atacLinuxRoot = '/mnt/a/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS'
$atacPython = '/home/micha/.local/share/atacseq-sept08-2026-venv/bin/python'

function Start-AtacStage {
    param([string]$Stage, [string]$Program, [string[]]$StageArguments)
    $stdout = Join-Path $atacRunRoot "logs/$Stage.stdout.log"
    $stderr = Join-Path $atacRunRoot "logs/$Stage.stderr.log"
    $process = Start-Process -FilePath $Program -ArgumentList $StageArguments -WorkingDirectory $atacRunRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $null = $process.Handle
    [PSCustomObject]@{Stage=$Stage; Process=$process; Stdout=$stdout; Stderr=$stderr}
}

function Wait-AtacStage {
    param($Job)
    $Job.Process.WaitForExit()
    $Job.Process.Refresh()
    if ($Job.Process.ExitCode -ne 0) {
        Get-Content -LiteralPath $Job.Stderr -Tail 30
        throw "Stage $($Job.Stage) failed with exit code $($Job.Process.ExitCode)"
    }
    Write-Output "Completed $($Job.Stage)"
}

try {
    if ($FinalizeOnly) {
        foreach ($required in @('counts/COMPLETE.json','logs/differential_session.txt',
            'logs/annotation_session.txt','logs/csaw_session.txt','logs/motifs_session.txt',
            'qc/matched_depth_peak_counts.tsv')) {
            if (-not (Test-Path -LiteralPath (Join-Path $atacRunRoot $required))) {
                throw "Cannot finalize: missing $required"
            }
        }
    } else {
    Write-Output 'Waiting for the completed consensus/count matrices.'
    while (-not (Test-Path -LiteralPath (Join-Path $atacRunRoot 'counts/COMPLETE.json'))) {
        Start-Sleep -Seconds 30
    }
    Write-Output 'Starting downstream ATAC analysis and matched-depth sensitivity.'
    $depthJob = Start-AtacStage -Stage '10_depth_sensitivity' -Program 'C:/Windows/System32/wsl.exe' -StageArguments @('-d','Ubuntu-24.04','--',$atacPython,'-u',"$atacLinuxRoot/code/10_depth_sensitivity.py")
    $differentialJob = Start-AtacStage -Stage '05_differential' -Program $atacRscript -StageArguments @('--vanilla',(Join-Path $atacRunRoot 'code/05_differential.R'))
    Wait-AtacStage $differentialJob
    $annotationJob = Start-AtacStage -Stage '06_annotation' -Program $atacRscript -StageArguments @('--vanilla',(Join-Path $atacRunRoot 'code/06_annotation.R'))
    Wait-AtacStage $annotationJob
    $csawJob = Start-AtacStage -Stage '07_csaw' -Program $atacRscript -StageArguments @('--vanilla',(Join-Path $atacRunRoot 'code/07_csaw.R'))
    $motifJob = Start-AtacStage -Stage '08_motifs' -Program $atacRscript -StageArguments @('--vanilla',(Join-Path $atacRunRoot 'code/08_motifs.R'))
    Wait-AtacStage $csawJob
    Wait-AtacStage $motifJob
    Wait-AtacStage $depthJob
    }
    $reportJob = Start-AtacStage -Stage '09_report' -Program 'C:/Windows/System32/wsl.exe' -StageArguments @('-d','Ubuntu-24.04','--',$atacPython,'-u',"$atacLinuxRoot/code/09_report.py")
    Wait-AtacStage $reportJob
    $verifyJob = Start-AtacStage -Stage '11_verify_outputs' -Program 'C:/Windows/System32/wsl.exe' -StageArguments @('-d','Ubuntu-24.04','--',$atacPython,'-u',"$atacLinuxRoot/code/11_verify_outputs.py")
    Wait-AtacStage $verifyJob
    Write-Output 'Report assembly and output verification completed. See the report for RNAseq/lpWGS integration status.'
} catch {
    $_ | Out-String | Set-Content -LiteralPath (Join-Path $atacRunRoot 'logs/downstream_error.txt')
    Write-Output $_
    exit 1
}
