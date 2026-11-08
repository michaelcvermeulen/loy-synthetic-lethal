$ErrorActionPreference = 'Stop'
$viewerPython = Join-Path $env:USERPROFILE 'anaconda3/python.exe'
if (-not (Test-Path -LiteralPath $viewerPython)) { $viewerPython = 'python' }
& $viewerPython (Join-Path $PSScriptRoot 'Start_interactive_report.py') --open
if ($LASTEXITCODE -ne 0) { throw 'Interactive report launch failed; see the error above.' }
