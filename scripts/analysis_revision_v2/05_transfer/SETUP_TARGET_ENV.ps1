param(
    [string]$Python = "py -3.12"
)

$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location $root

if (-not (Test-Path -LiteralPath ".\analysis\.venv\Scripts\python.exe")) {
    Invoke-Expression "$Python -m venv .\analysis\.venv"
}

$venvPython = ".\analysis\.venv\Scripts\python.exe"
& $venvPython -m pip install --upgrade pip
& $venvPython -m pip install --requirement .\revision_v2\05_transfer\requirements-gate03.txt

$env:GSE206528_SCANPY_SNAPSHOT_DIR = (Resolve-Path ".\tools\scanpy_snapshot").Path
& $venvPython .\revision_v2\05_transfer\VERIFY_TARGET_READY.py --write-report
if ($LASTEXITCODE -ne 0) {
    throw "Target readiness check failed. Do not start the replay."
}

Write-Host "Environment installed and target preflight passed."
