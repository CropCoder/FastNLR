#!/usr/bin/env bash
#
# 下载并安装 FastNLR 最新版预编译二进制到 ~/.local/bin。
#
# 用法：bash install_fastnlr.sh [版本号，如 v2.0.0]
#
# Author: Jiwen Zhao (https://github.com/CropCoder)

set -euo pipefail

version="${1:-v2.0.0}"
base_url="https://github.com/CropCoder/FastNLR/releases/download/${version}"

case "$(uname -s)-$(uname -m)" in
  Linux-x86_64) target="x86_64-unknown-linux-gnu" ;;
  Darwin-x86_64) target="x86_64-apple-darwin" ;;
  Darwin-arm64) target="aarch64-apple-darwin" ;;
  *)
    echo "unsupported platform: $(uname -s)-$(uname -m)" >&2
    exit 1
    ;;
esac

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

archive="fastnlr-${target}.tar.gz"
curl -fL "${base_url}/${archive}" -o "${tmpdir}/${archive}"
curl -fL "${base_url}/${archive}.sha256" -o "${tmpdir}/${archive}.sha256"

(
  cd "${tmpdir}"
  sha256sum -c "${archive}.sha256" >/dev/null
  tar -xzf "${archive}"
)

install_dir="${CODEX_INSTALL_DIR:-$HOME/.local/bin}"
mkdir -p "${install_dir}"
install -m 0755 "${tmpdir}/fastnlr" "${install_dir}/fastnlr"
echo "installed fastnlr ${version} to ${install_dir}/fastnlr"
