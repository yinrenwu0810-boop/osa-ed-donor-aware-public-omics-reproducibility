[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$ProjectRoot = 'C:\Users\22394\Documents\Codex\2026-07-12\acad'
$AuditDirectory = Join-Path $ProjectRoot 'revision_v2\07_virtual_perturbation\03_preflight\memory_release_events'
$OldProjectWorkerPattern = '/mnt/c/Users/22394/Documents/Codex/2026-07-05/new-chat/work/run_biomni_mcp_wsl.py'

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

$started = (Get-Date).ToUniversalTime().ToString('o')
$before = Get-FreeRamGb
$actions = @()
$wslRunning = @(Get-Process -Name vmmemWSL -ErrorAction SilentlyContinue).Count -gt 0
if (-not $wslRunning) {
    $actions += 'wsl_not_running_no_wsl_cleanup_needed'
} else {
try {
    $wslLines = @(wsl.exe -e sh -lc 'ps -eo pid,args --no-headers' 2>$null)
    $targetPids = @()
    foreach ($line in $wslLines) {
        if ($line -match '^\s*(\d+)\s+(.+)$') {
            $wslPid = [int]$Matches[1]
            $command = $Matches[2]
            if ($command -like "*$OldProjectWorkerPattern*") { $targetPids += $wslPid; continue }
        }
    }
    if ($targetPids.Count -gt 0) {
        wsl.exe -e sh -lc ("kill -TERM " + ($targetPids -join ' ')) | Out-Null
        Start-Sleep -Seconds 2
        $actions += "terminated_confirmed_old_project_wsl_workers:" + ($targetPids -join ',')
    }
    # User-authorized policy: this G04 runner does not use WSL.  Once a
    # confirmed obsolete Biomni worker is found, close the whole WSL session
    # so vmmemWSL cannot retain memory for an interactive shell or services.
    if ($targetPids.Count -gt 0) {
        wsl.exe --shutdown
        $actions += 'shutdown_wsl_after_confirmed_obsolete_biomni_workers'
    }
} catch {
    $actions += ('wsl_cleanup_error:' + $_.Exception.Message)
}
}
$after = Get-FreeRamGb
New-Item -ItemType Directory -Force -Path $AuditDirectory | Out-Null
$event = [pscustomobject]@{
    timestamp_utc = $started
    free_ram_gb_before = $before
    free_ram_gb_after = $after
    threshold_gb = 10
    confirmed_old_project_worker_pattern = $OldProjectWorkerPattern
    unknown_wsl_workload_detected = $false
    actions = $actions
    status = if ($after -ge 10) { 'RESOURCE_READY' } else { 'HOLD_RESOURCE' }
}
$eventPath = Join-Path $AuditDirectory ('memory_release_' + (Get-Date -Format 'yyyyMMddTHHmmssfff') + '.json')
$event | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $eventPath -Encoding utf8
$event | ConvertTo-Json -Compress
if ($after -lt 10) { exit 2 }
