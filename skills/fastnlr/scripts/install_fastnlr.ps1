# 下载并安装 FastNLR 最新版预编译二进制到当前用户目录。
#
# 用法：powershell -ExecutionPolicy Bypass -File install_fastnlr.ps1 [-Version v1.7.1]
#
# Author: Jiwen Zhao (https://github.com/CropCoder)

param([string]$Version = "v1.7.1")

$ErrorActionPreference = "Stop"
$base = "https://github.com/CropCoder/FastNLR/releases/download/$Version"
$target = "x86_64-pc-windows-msvc"
$zip = "fastnlr-$target.zip"
$tmp = Join-Path $env:TEMP "fastnlr-$Version.zip"

Invoke-WebRequest -Uri "$base/$zip" -OutFile $tmp
$dir = Join-Path $env:LOCALAPPDATA "FastNLR"
New-Item -ItemType Directory -Force -Path $dir | Out-Null
Expand-Archive -Path $tmp -DestinationPath $dir -Force

Write-Host "installed fastnlr $Version to $dir"
