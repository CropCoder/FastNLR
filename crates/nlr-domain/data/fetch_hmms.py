#!/usr/bin/env python3
"""Download the curated Pfam HMMs used by FastNLR protein mode.

用法:
    python3 fetch_hmms.py [输出目录]      # 默认输出到脚本所在目录

说明:
    - 数据来源: InterPro REST API 的 Pfam HMM 注解接口，返回 gzip 压缩的 HMMER3 ASCII 文件；
    - 脚本会解压为 `.hmm` 文本，文件名即 Pfam accession；
    - 版本/阈值会随 Pfam 更新，重新下载后请运行 `cargo test -p nlr-domain`
      以确认 golden-score 断言（±0.5 bit）是否需要同步更新。
"""

import gzip
import io
import os
import sys
import urllib.request

MODELS = ["PF00931", "PF01582", "PF18052", "PF05659"]
URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}/?annotation=hmm"


def fetch(acc: str) -> bytes:
    url = URL.format(acc=acc)
    with urllib.request.urlopen(url, timeout=60) as resp:
        raw = resp.read()
    try:
        return gzip.decompress(raw)
    except OSError:
        return raw  # 已经是未压缩文本时直接使用


def main() -> int:
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)
    for acc in MODELS:
        text = fetch(acc)
        header = text.split(b"\n", 1)[0]
        if not header.startswith(b"HMMER3"):
            print(f"warning: {acc} 的响应不像 HMMER3 文件: {header[:40]!r}", file=sys.stderr)
        path = os.path.join(outdir, f"{acc}.hmm")
        with open(path, "wb") as fh:
            fh.write(text)
        print(f"{path}  ({len(text)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
