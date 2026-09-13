[CmdletBinding()]
param([int]$StartBatchOrder = 1)
$ErrorActionPreference = 'Stop'
$ProjectRoot = 'C:\Users\22394\Documents\Codex\2026-07-12\acad'
$VPRoot = Join-Path $ProjectRoot 'revision_v2\07_virtual_perturbation'
$Protocol = Join-Path $VPRoot '01_protocol\control_execution_v2'
$Matrix = Join-Path $Protocol 'control_batch_matrix.tsv'
$Freeze = Join-Path $Protocol 'VP_G05_control_execution_freeze_v2.json'
$Runner = Join-Path $VPRoot 'scripts\05d_run_control_batch_shared_wt.R'
$RScript = 'D:\R\R-4.6.1\bin\Rscript.exe'
$Active = Join-Path $VPRoot '05_runs\active_task_VP_G05_control_batch.json'
$Cleanup = Join-Path $VPRoot 'scripts\04_release_unrelated_memory.ps1'
if (-not (Test-Path -LiteralPath $Freeze)) { throw 'Control execution is not frozen.' }
if (Test-Path -LiteralPath $Active) { throw "An active control batch is already recorded: $Active" }
$freezeHash = (Get-FileHash -LiteralPath $Freeze -Algorithm SHA256).Hash.ToLowerInvariant()
function Get-FreeRamGb {
    try { $os = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop; return [math]::Round($os.FreePhysicalMemory / 1MB, 2) }
    catch { $bytes = @(& python -c 'import psutil; print(psutil.virtual_memory().available)' 2>$null | Select-Object -Last 1); if ($LASTEXITCODE -ne 0 -or $bytes[0] -notmatch '^\d+$') { throw 'Unable to determine available RAM.' }; return [math]::Round(([double]$bytes[0] / 1GB), 2) }
}
$batches = Import-Csv -LiteralPath $Matrix -Delimiter "`t" | Sort-Object {[int]$_.batch_order} | Where-Object {[int]$_.batch_order -ge $StartBatchOrder}
foreach ($batch in $batches) {
    $outputBatch = Join-Path $VPRoot ($batch.output_batch -replace '/', '\')
    if (Test-Path -LiteralPath (Join-Path $outputBatch 'batch_complete.json')) { continue }
    $free = Get-FreeRamGb
    if ($free -lt 10) { & $Cleanup | Out-Host; $free = Get-FreeRamGb }
    if ($free -lt 10) { throw "HOLD_RESOURCE: FreeRAM_GB=$free below 10 GB before batch $($batch.batch_order)." }
    $controls = Join-Path $Protocol ($batch.controls_file -replace '/', '\')
    $matrixDir = Join-Path $VPRoot ($batch.matrix_directory -replace '/', '\')
    $benchmark = Join-Path $VPRoot ($batch.benchmark_file -replace '/', '\')
    $activeRecord = [ordered]@{gate='VP-G05'; stage='matched_control_KO_batch'; status='RUNNING'; batch_order=[int]$batch.batch_order; condition=$batch.condition; donor=$batch.donor; seed=[int]$batch.seed; controls=[int]$batch.control_count; output_batch=$outputBatch; started_at_utc=(Get-Date).ToUniversalTime().ToString('o'); launcher_pid=$PID}
    $activeRecord | ConvertTo-Json | Set-Content -LiteralPath $Active -Encoding utf8
    & $RScript $Runner --matrix-dir $matrixDir --controls-file $controls --benchmark-file $benchmark --output-batch $outputBatch --condition $batch.condition --donor $batch.donor --benchmark-gene $batch.benchmark_gene --seed $batch.seed --n-cores 4 --freeze-sha256 $freezeHash
    $exitCode = $LASTEXITCODE
    Remove-Item -LiteralPath $Active -Force -ErrorAction SilentlyContinue
    if ($exitCode -ne 0) { throw "Control batch $($batch.batch_order) failed or held with exit code $exitCode." }
}
