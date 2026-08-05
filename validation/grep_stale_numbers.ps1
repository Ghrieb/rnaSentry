# Literal repo-wide sweep for superseded concordance numbers and the old
# "poor generalization" narrative. Scans .R (code + @examples), .Rmd (vignette
# and report template), .md, and .Rd (man pages), excluding build artifacts
# (*.Rcheck/ and chk*/ directories). Hits are only permitted under validation/
# and the case-study narrative vignettes (the intentional before/after
# correction narrative and the Tier-1 naive-vs-guarded quotes in the case
# studies); any hit elsewhere fails the gate and exits non-zero.
#
# Usage (from the rnaSentry package root):
#   powershell -ExecutionPolicy Bypass -File validation/grep_stale_numbers.ps1

param(
  [string]$Root = (Split-Path -Parent $PSScriptRoot)
)

$patterns = '0\.217|0\.424|0\.510|0\.160|poor generalization|winner.{0,3}s curse|downward selection|0\.22 collapse'
$exts = '.R', '.Rmd', '.md', '.Rd'

# Allowed correction narrative: the validation dossier plus the case-study
# vignettes that deliberately quote the superseded numbers as the naive
# counterfactual (case-study-brca and the case-study-impact landing page).
$allowed = @(
  'validation\*',
  'vignettes\case-study-brca.Rmd',
  'vignettes\case-study-impact.Rmd'
)

$files = Get-ChildItem -Path $Root -Recurse -File -Force |
  Where-Object {
    $exts -contains $_.Extension -and
    $_.FullName -notmatch '\\.git\\' -and
    $_.FullName -notmatch '\.Rcheck\\' -and
    $_.FullName -notmatch '\\chk[0-9]+\\'
  }

$hits = @()
foreach ($f in $files) {
  $m = Select-String -LiteralPath $f.FullName -Pattern $patterns
  foreach ($x in $m) {
    $rel = $f.FullName.Substring($Root.Length).TrimStart('\', '/')
    $hits += [pscustomobject]@{ File = $rel; Line = $x.LineNumber; Text = $x.Line.Trim() }
  }
}

Write-Host "== Stale-number sweep =="
Write-Host "Pattern(s): $patterns"
Write-Host "Files scanned: $($files.Count)  (.R/.Rmd/.md/.Rd)"
if ($hits.Count -eq 0) {
  Write-Host "No hits."
} else {
  $hits | Sort-Object File, Line |
    ForEach-Object { "{0}:{1}: {2}" -f $_.File, $_.Line, $_.Text }
  Write-Host ""
  Write-Host "Hits in allowed correction narrative: $($hits.Where({ $f = $_.File; $allowed | Where-Object { $f -like $_ } }).Count)"
}

$unexpected = $hits | Where-Object {
  $f = $_.File
  -not ($allowed | Where-Object { $f -like $_ })
}
if ($unexpected.Count -gt 0) {
  Write-Host "FAIL: $($unexpected.Count) stale reference(s) outside the allowed narrative." -ForegroundColor Red
  exit 1
}
Write-Host "PASS: no stale references outside the allowed narrative." -ForegroundColor Green
exit 0
