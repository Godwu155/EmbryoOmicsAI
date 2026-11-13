param(
    [Parameter(Mandatory = $true)][string]$Exchange,
    [string]$Environment = 'results/cellxgene_env_20261113',
    [ValidateRange(1024, 65535)][int]$Port = 5057
)
$ErrorActionPreference = 'Stop'
$taskRepo = Split-Path $PSScriptRoot -Parent
$taskEnvironment = if ([IO.Path]::IsPathRooted($Environment)) { $Environment } else { Join-Path $taskRepo $Environment }
$taskExecutable = Join-Path $taskEnvironment 'Scripts/cellxgene.exe'
$taskExchange = (Resolve-Path -LiteralPath $Exchange).Path
$taskViewer = Join-Path $taskExchange 'viewer.h5ad'
$taskLabels = Join-Path $taskExchange 'cell_annotations.csv'
$taskManifest = Get-Content -LiteralPath (Join-Path $taskExchange 'exchange.json') -Raw | ConvertFrom-Json
if (-not (Test-Path -LiteralPath $taskExecutable)) { throw 'Install the separate CELLxGENE environment described in docs/CELLXGENE_ACCEPTANCE.md' }
if ((Get-FileHash -LiteralPath $taskViewer -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskManifest.viewer_sha256) { throw 'Viewer file changed since export; export a fresh browsing copy' }
if (-not (Test-Path -LiteralPath $taskLabels)) { throw 'Missing cell_annotations.csv' }
& $taskExecutable launch $taskViewer --host 127.0.0.1 --port $Port --var-names viewer_gene --annotations-file $taskLabels --disable-diffexp --disable-gene-sets-save --max-category-items 500 --title 'EmbryoOmics exploratory review'
exit $LASTEXITCODE
