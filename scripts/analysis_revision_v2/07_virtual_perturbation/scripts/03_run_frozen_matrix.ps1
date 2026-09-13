[CmdletBinding()]
param(
    [int]$StartRunOrder = 1,
    [int]$MaxRuns = 45
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = 'C:\Users\22394\Documents\Codex\2026-07-12\acad'
$VPRoot = Join-Path $ProjectRoot 'revision_v2\07_virtual_perturbation'
$RScript = 'D:\R\R-4.6.1\bin\Rscript.exe'
$RLib = Join-Path $VPRoot '02_env\Rlib'
$RunMatrix = Join-Path $VPRoot '01_protocol\run_matrix.tsv'
$Runner = Join-Path $VPRoot 'scripts\02_run_sctenifoldknk.R'
$RunsRoot = Join-Path $VPRoot '05_runs\main'
$Registry = Join-Path $VPRoot '05_runs\run_registry.tsv'
$ActiveTask = Join-Path $VPRoot '05_runs\active_task.json'
$MemoryCleanup = Join-Path $VPRoot 'scripts\04_release_unrelated_memory.ps1'
$PauseEvents = Join-Path $VPRoot '05_runs\pause_events'

if (-not (Test-Path -LiteralPath (Join-Path $VPRoot 'validation\VP_G03_pilot_PASS.json'))) { throw 'VP-G03 is not PASS.' }
if (Test-Path -LiteralPath $ActiveTask) { throw "An active G04 task is already recorded: $ActiveTask" }
if (-not (Test-Path -LiteralPath $RScript) -or -not (Test-Path -LiteralPath $RLib)) { throw 'Frozen runner or task-specific R library is unavailable.' }
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Registry), $RunsRoot | Out-Null
$env:R_LIBS_USER = $RLib

function Get-FreeRamGb {
    try {
        $os = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
        return [math]::Round($os.FreePhysicalMemory / 1MB, 2)
    } catch {
        # The managed desktop sandbox can deny WMI queries.  psutil reads the
        # same OS memory counter without broadening process-management scope.
        $availableBytes = @(& python -c 'import psutil; print(psutil.virtual_memory().available)' 2>$null | Select-Object -Last 1)
        if ($LASTEXITCODE -ne 0 -or $availableBytes.Count -ne 1 -or $availableBytes[0] -notmatch '^\d+$') {
            throw 'Unable to determine available RAM through WMI or psutil.'
        }
        return [math]::Round(([double]$availableBytes[0] / 1GB), 2)
    }
}

function Add-FinalRegistryRow {
    param([pscustomobject]$Task, [string]$Status, [int]$ExitCode, [string]$AttemptPath, [string]$Note, [string]$StartedAt)
    $row = [pscustomobject]@{
        run_order = [int]$Task.run_order
        run_id = $Task.run_id
        gene = $Task.gene
        condition = $Task.condition
        donor = $Task.donor
        seed = [int]$Task.seed
        started_at_utc = $StartedAt
        finished_at_utc = (Get-Date).ToUniversalTime().ToString('o')
        status = $Status
        exit_code = $ExitCode
        output_attempt = $AttemptPath
        note = $Note
    }
    $row | Export-Csv -LiteralPath $Registry -NoTypeInformation -Append -Encoding utf8
}

$tasks = Import-Csv -LiteralPath $RunMatrix -Delimiter "`t" | Sort-Object { [int]$_.run_order } | Where-Object { [int]$_.run_order -ge $StartRunOrder } | Select-Object -First $MaxRuns
foreach ($task in $tasks) {
    $completed = @(Import-Csv -LiteralPath $Registry | Where-Object { [int]$_.run_order -eq [int]$task.run_order -and $_.status -eq 'PASS_TECHNICAL' })
    if ($completed.Count -gt 0) { continue }
    $freeRamGb = Get-FreeRamGb
    if ($freeRamGb -lt 10) {
        & $MemoryCleanup
        $freeRamGb = Get-FreeRamGb
    }
    if ($freeRamGb -lt 10) {
        Add-FinalRegistryRow -Task $task -Status 'HOLD_RESOURCE' -ExitCode 0 -AttemptPath '' -Note "FreeRAM_GB=$freeRamGb below 10 GB before task launch; no model run started." -StartedAt ''
        break
    }
    $matrixDir = Join-Path $VPRoot ($task.matrix_directory -replace '/', '\\')
    $targetQc = @(Import-Csv -LiteralPath (Join-Path $matrixDir 'target_detection.tsv') -Delimiter "`t" | Where-Object { $_.gene -eq $task.gene })
    if ($targetQc.Count -ne 1 -or [string]$targetQc[0].status -ne 'PASS') {
        Add-FinalRegistryRow -Task $task -Status 'HOLD_TARGET_DETECTION' -ExitCode 0 -AttemptPath '' -Note 'Target detection QC is not PASS; no model run started.' -StartedAt ''
        break
    }
    $attemptPath = Join-Path $VPRoot ($task.output_attempt_directory -replace '/', '\\')
    if (Test-Path -LiteralPath $attemptPath) {
        $pauseEvidence = @(Get-ChildItem -LiteralPath $PauseEvents -Filter ($task.run_id + '_PAUSED_*.json') -File -ErrorAction SilentlyContinue)
        if ($pauseEvidence.Count -eq 0) {
            Add-FinalRegistryRow -Task $task -Status 'FAIL_RETAINED' -ExitCode 0 -AttemptPath $attemptPath -Note 'Attempt path already exists before launch; refusing overwrite.' -StartedAt ''
            break
        }
        $attemptParent = Split-Path -Parent $attemptPath
        $attemptNumber = 2
        do {
            $attemptPath = Join-Path $attemptParent ('attempt_{0:D2}' -f $attemptNumber)
            $attemptNumber++
        } while (Test-Path -LiteralPath $attemptPath)
    }
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $attemptPath) | Out-Null
    $startedAt = (Get-Date).ToUniversalTime().ToString('o')
    @{ run_order = [int]$task.run_order; run_id = $task.run_id; gene = $task.gene; condition = $task.condition; donor = $task.donor; seed = [int]$task.seed; matrix_dir = $matrixDir; output_attempt = $attemptPath; started_at_utc = $startedAt; launcher_pid = $PID } | ConvertTo-Json | Set-Content -LiteralPath $ActiveTask -Encoding utf8
    & $RScript $Runner --matrix-dir $matrixDir --target $task.gene --seed $task.seed --output-attempt $attemptPath --n-cores 4
    $exitCode = $LASTEXITCODE
    Remove-Item -LiteralPath $ActiveTask -Force -ErrorAction SilentlyContinue
    $requiredFiles = @('differential_regulation.tsv','run_parameters.json','session_info.txt','stdout.log','stderr.log','resource_usage.tsv','artifact_manifest.sha256.tsv')
    $technicalPass = $exitCode -eq 0 -and (Test-Path -LiteralPath $attemptPath) -and (($requiredFiles | Where-Object { -not (Test-Path -LiteralPath (Join-Path $attemptPath $_)) }).Count -eq 0)
    if ($technicalPass) {
        Add-FinalRegistryRow -Task $task -Status 'PASS_TECHNICAL' -ExitCode $exitCode -AttemptPath $attemptPath -Note 'Frozen task completed; no biological interpretation performed.' -StartedAt $startedAt
        continue
    }
    $failurePath = "$attemptPath`_FAIL"
    $recordPath = if (Test-Path -LiteralPath $failurePath) { $failurePath } elseif (Test-Path -LiteralPath $attemptPath) { $attemptPath } else { '' }
    Add-FinalRegistryRow -Task $task -Status 'FAIL_RETAINED' -ExitCode $exitCode -AttemptPath $recordPath -Note 'Runner returned nonzero or required technical artifacts were absent.' -StartedAt $startedAt
    break
}
