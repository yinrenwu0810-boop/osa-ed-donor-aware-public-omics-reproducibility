$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location $root
$venvPython = ".\analysis\.venv\Scripts\python.exe"
if (-not (Test-Path -LiteralPath $venvPython)) {
    throw "Run SETUP_TARGET_ENV.ps1 first."
}

$env:GSE206528_SCANPY_SNAPSHOT_DIR = (Resolve-Path ".\tools\scanpy_snapshot").Path
& $venvPython .\revision_v2\05_transfer\VERIFY_TARGET_READY.py --write-report
if ($LASTEXITCODE -ne 0) {
    throw "Target preflight did not pass; the replay was not started."
}

& $venvPython .\revision_v2\scripts\run_GSE206528_repro_v2.py
$runnerExit = $LASTEXITCODE
$attempt = Get-ChildItem -LiteralPath ".\revision_v2\01_work\GSE206528_replay" -Directory |
    Sort-Object Name | Select-Object -Last 1

if ($null -eq $attempt) {
    throw "The runner did not create an attempt directory."
}

& $venvPython .\revision_v2\scripts\check_GSE206528_replay_v2.py $attempt.FullName
$checkerExit = $LASTEXITCODE

if ($runnerExit -ne 0 -or $checkerExit -ne 0) {
    exit 2
}
