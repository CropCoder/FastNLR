#!/usr/bin/env bash
#
# 打包 FastNLR Web 服务为 TBtools 可拉起的本地 Web 插件。
#
# 用法：bash packaging/Build-TBtools-Plugin-V1.sh <target> [version]
#   target  linux | macos | windows
#   version 默认 v2.0.0
#
# 说明：
#   - 需在目标操作系统上运行（macOS/Windows 分别打包）。
#   - 产物为 dist/fastnlr-server-tbtools-<target>.plugin.zip。
#   - 正式 .plugin 请用 TBtools 的 Plugin Creator 包装本目录的二进制与 launch 脚本。
#
# Author: Jiwen Zhao (https://github.com/CropCoder)

set -euo pipefail

target="${1:?用法: Build-TBtools-Plugin-V1.sh <linux|macos|windows> [version]}"
version="${2:-v2.0.0}"

case "$target" in
  linux)   rust_target="x86_64-unknown-linux-gnu"; ext="" ;;
  macos)   rust_target="x86_64-apple-darwin"; ext="" ;;
  windows) rust_target="x86_64-pc-windows-msvc"; ext=".exe" ;;
  *) echo "unknown target: $target" >&2; exit 1 ;;
esac

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
out_dir="${repo_root}/dist/tbtools-${target}"
mkdir -p "$out_dir"

cargo build --release -p nlr-server --target "$rust_target"

bin_name="fastnlr-server${ext}"
cp "${repo_root}/target/${rust_target}/release/${bin_name}" "${out_dir}/${bin_name}"

cat > "${out_dir}/plugin.json" <<EOF
{
  "name": "FastNLR",
  "version": "${version}",
  "description": "Fast, self-contained NLR immune-receptor locus annotation for plant genomes",
  "type": "web",
  "launch": "${bin_name}",
  "url": "http://127.0.0.1:8787"
}
EOF

cat > "${out_dir}/launch.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
port="${FASTNLR_PORT:-8787}"
"$(dirname "$0")/fastnlr-server" --host 127.0.0.1 --port "$port" &
sleep 1
open "http://127.0.0.1:${port}" 2>/dev/null || xdg-open "http://127.0.0.1:${port}" 2>/dev/null || true
wait
EOF

cat > "${out_dir}/launch.bat" <<'EOF'
@echo off
set PORT=8787
start /b fastnlr-server.exe --host 127.0.0.1 --port %PORT%
start http://127.0.0.1:%PORT%
EOF

chmod +x "${out_dir}/launch.sh"

zip_name="fastnlr-server-tbtools-${target}.plugin.zip"
(cd "${out_dir}" && zip -qr "../${zip_name}" .)
echo "built: ${repo_root}/dist/${zip_name}"
