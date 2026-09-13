$ErrorActionPreference = 'Stop'

$inputPath = Join-Path $PSScriptRoot 'TYMS_EFNB2_LRRC17_虚拟扰动实施工作手册_v1.docx'
$renderDir = Join-Path $PSScriptRoot '_qa_render_word_v2'
$outputPath = Join-Path $renderDir 'TYMS_EFNB2_LRRC17_虚拟扰动实施工作手册_v1.pdf'

New-Item -ItemType Directory -Path $renderDir -Force | Out-Null

$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    $document = $word.Documents.Open([System.IO.Path]::GetFullPath($inputPath), $false, $true)
    try {
        $document.ExportAsFixedFormat([System.IO.Path]::GetFullPath($outputPath), 17)
    }
    finally {
        $document.Close($false)
    }
    Write-Output $outputPath
}
finally {
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}
