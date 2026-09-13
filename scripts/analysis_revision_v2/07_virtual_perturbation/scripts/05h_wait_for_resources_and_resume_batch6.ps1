[CmdletBinding()]
param(
    [double]$MinimumFreeRamGb = 10,
    [int]$PollSeconds = 60,
    [int]$MaximumWaitHours = 24
)
$ErrorActionPreference = 'Stop'
$ProjectRoot = 'C:\Users\22394\Documents\Codex\2026-07-12\acad'
$VPRoot = Join-Path $ProjectRoot 'revision_v2\07_virtual_perturbation'
$Launcher = Join-Path $VPRoot 'scripts\05f2_run_control_batches_v2_compat.ps1'
$ActiveBatch = Join-Path $VPRoot '05_runs\active_task_VP_G05_control_batch.json'
$CompletedBatch6 = Join-Path $VPRoot '05_runs\control_batches_v2\organic_ED_nonDM\non-DM_3\2026082701\attempt_01\batch_complete.json'
$WaitLog = Join-Path $VPRoot '05_runs\VP_G05_batch6_resource_wait.tsv'

function Get-FreeRamGb {
    try {
        $os = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
        return [math]::Round($os.FreePhysicalMemory / 1MB, 2)
    } catch {
        $bytes = @(& python -c 'import psutil; print(psutil.virtual_memory().available)' 2>$null | Select-Object -Last 1)
        if ($LASTEXITCODE -ne 0 -or $bytes.Count -ne 1 -or $bytes[0] -notmatch '^\d+$') { throw 'Unable to determine available RAM.' }
        return [math]::Round(([double]$bytes[0] / 1GB), 2)
    }
}
function Add-WaitRecord([string]$Status, [double]$FreeRamGb) {
    $record = [pscustomobject]@{ timestamp_utc=(Get-Date).ToUniversalTime().ToString('o'); status=$Status; free_ram_gb=$FreeRamGb; threshold_gb=$MinimumFreeRamGb; monitor_pid=$PID }
    $record | Export-Csv -LiteralPath $WaitLog -Delimiter "`t" -NoTypeInformation -Append -Encoding utf8
}

if (Test-Path -LiteralPath $CompletedBatch6) { Add-WaitRecord 'ALREADY_COMPLETE' (Get-FreeRamGb); exit 0 }
$deadline = (Get-Date).ToUniversalTime().AddHours($MaximumWaitHours)
while ((Get-Date).ToUniversalTime() -lt $deadline) {
    if (Test-Path -LiteralPath $CompletedBatch6) { Add-WaitRecord 'ALREADY_COMPLETE' (Get-FreeRamGb); exit 0 }
    if ((Test-Path -LiteralPath $ActiveBatch) -or (Get-Process Rscript,Rterm -ErrorAction SilentlyContinue)) { Add-WaitRecord 'EXIT_OTHER_RUN_DETECTED' (Get-FreeRamGb); exit 3 }
    $free = Get-FreeRamGb
    if ($free -ge $MinimumFreeRamGb) {
        Add-WaitRecord 'RESOURCE_READY_STARTING_BATCH6' $free
        & $Launcher -StartBatchOrder 6
        exit $LASTEXITCODE
    }
    Add-WaitRecord 'WAITING_RESOURCE' $free
    Start-Sleep -Seconds $PollSeconds
}
Add-WaitRecord 'HOLD_RESOURCE_WAIT_TIMEOUT' (Get-FreeRamGb)
exit 2
