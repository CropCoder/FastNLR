#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把 MEME 挖出的 RNL（helper NLR：ADR1/NRG1/RPW8 类）motif 转成 FastNLR 的 mot.txt / store.txt 条目。

FastNLR 的 motif 机制（读源码 nlr-config / nlr-scan 得到）：
  · mot.txt : `motif_<id>@<pos0based>@<LETTER> <int>`，整数得分，字母表 A–Z（索引 = ascii-65）；
               现有 20 个 motif 的取值都在 0–100。
  · store.txt: `motif_<id>@<score> <p>`，**右尾 CDF**（p = P(随机窗口得分 ≥ score)），单调递减；
               扫描时把窗口内各位置得分求和，用这个和当索引查 p；p < 1e-5 才算命中。
  · 索引上限 = 词宽 × 100（现有 motif 的 CDF 长度正好是 width*100+1）。
  · 另外 CLI 还有一道预过滤 `has_signature`：一个 fragment 的 motif 串里必须能匹配到 SIGNATURES 里的
    连续子串，否则该 fragment 的命中被整体丢弃 —— 所以新增 motif 也要进 SIGNATURES。

本脚本做三件事：
  1) 解析 MEME 的 letter-probability matrix → 整数 PWM（p×100，保持与现有一致的 0–100 量纲）；
  2) 用**精确动态规划**（逐位置卷积 + 大豆蛋白组的氨基酸背景频率）算出窗口得分分布，
     写成右尾 CDF —— 不用随机采样，尾部精度到 1e-12 级别；
  3) 输出可直接追加进 mot.txt / store.txt 的文本，并报告 p=1e-4 / 1e-5 对应的得分阈值与
     在蛋白组里的预期命中数（用来判断新 motif 的特异性）。

用法：
  python3 make_rnl_motifs.py --meme /tmp/meme_rnl/meme.txt --proteome work_t2t/proteome.faa \
      --start-id 21 --out-prefix /tmp/rnl_motif
"""
import argparse
import collections
import os
import sys

AA20 = 'ACDEFGHIKLMNPQRSTVWY'
# mot.txt 的字母表是 A–Z（索引 = ascii-65）：除 20 个标准氨基酸外还有 B/Z/J/X/U/O
LETTERS = [chr(c) for c in range(ord('A'), ord('Z') + 1)]


def read_fasta(path):
    seqs, name, buf = [], None, []
    op = open
    if path.endswith('.gz'):
        import gzip
        op = gzip.open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    seqs.append(''.join(buf))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        seqs.append(''.join(buf))
    return seqs


def parse_meme(path):
    """→ [(name, width, [[p_aa per position] for 20 aa])]

    MEME 文本格式：`MOTIF <consensus> MEME-<n>` 行开始，`letter-probability matrix: alength= 20 w= W ...`
    之后紧跟 W 行、每行 20 个浮点数（字母序 A C D E F G H I K L M N P Q R S T V W Y）。
    """
    motifs = []
    name = None
    mode = None
    width = 0
    matrix = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('MOTIF '):
                if name and matrix:
                    motifs.append((name, width, matrix))
                name = line.split()[1]
                mode, width, matrix = None, 0, []
                continue
            if 'letter-probability matrix' in line:
                mode = 'matrix'
                width = int(line.split('w=')[1].split()[0])
                matrix = []
                continue
            if mode == 'matrix' and name:
                parts = line.split()
                if len(parts) == 20 and _isfloat(parts[0]):
                    matrix.append([float(x) for x in parts])
                    if len(matrix) == width:
                        mode = None
                elif matrix:
                    mode = None
    if name and matrix:
        motifs.append((name, width, matrix))
    return motifs


def _isfloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def bg_freq(seqs, pseudocount=0.0):
    c = collections.Counter()
    n = 0
    for s in seqs:
        for a in s:
            if a in AA20:
                c[a] += 1
                n += 1
    freqs = {a: (c[a] + pseudocount) / (n + 20 * pseudocount) for a in AA20}
    return freqs, n


def score_matrix(matrix, scale=100.0):
    """MEME 概率矩阵 → 整数 PWM（0..scale，与现有 0–100 量纲一致）"""
    return [[int(round(p * scale)) for p in row] for row in matrix]


def cdf_exact(smatrix, freqs):
    """精确窗口得分分布：逐位置卷积（每位置得分 = score[pos][aa]，权重 = 背景频率）

    → (cdf 列表, pmf 列表)：cdf[s] = P(窗口得分 ≥ s)
    """
    import numpy as np
    maxs = 100 * len(smatrix)
    pmf = np.zeros(maxs + 1)
    pmf[0] = 1.0
    for row in smatrix:
        step = np.zeros(maxs + 1)
        for a, sc in zip(AA20, row):
            step[sc] += freqs[a]
        pmf = np.convolve(pmf, step)[:maxs + 1]
    tail = np.cumsum(pmf[::-1])[::-1]        # 右尾累积
    return tail.tolist(), pmf.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--meme', required=True)
    ap.add_argument('--proteome', required=True)
    ap.add_argument('--start-id', type=int, default=21)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--scale', type=float, default=100.0)
    ap.add_argument('--keep', default='', help='只保留 MEME 里第几个 motif（1-based，逗号分隔），'
                                             '如 1,2,3；留空=全要')
    ap.add_argument('--target-exp', type=float, default=None,
                    help='收紧标定：把 CDF 右移，使"全蛋白组期望命中窗口数（p<1e-5）"不超过该值。'
                         'TIR 家族 motif 命中面很宽（默认每种 motif 约 370 个期望窗口），'
                         '实测不收紧会让全基因组位点数从 690 涨到 1168、且多为<1kb 的短位点；'
                         '收紧后只让最强的命中播种（同时保留单子叶 TIR-only 的检出）')
    a = ap.parse_args()

    motifs = parse_meme(a.meme)
    if a.keep:
        keep = {int(x) for x in a.keep.split(',') if x.strip()}
        motifs = [m for i, m in enumerate(motifs, 1) if i in keep]
    print('MEME motif 数：%d%s' % (len(motifs), ('（保留 %s）' % a.keep) if a.keep else ''), file=sys.stderr)
    seqs = read_fasta(a.proteome)
    freqs, nres = bg_freq(seqs)
    print('背景：%d 条蛋白 / %.3fM 残基' % (len(seqs), nres / 1e6), file=sys.stderr)

    mot_out, store_out, report = [], [], []
    for k, (name, width, matrix) in enumerate(motifs):
        mid = a.start_id + k
        smatrix = score_matrix(matrix, a.scale)
        tail, pmf = cdf_exact(smatrix, freqs)
        # 阈值
        def thr(p):
            for s, v in enumerate(tail):
                if v < p:
                    return s
            return len(tail)
        t4, t5 = thr(1e-4), thr(1e-5)
        # 逐位置写 PWM：20 个标准氨基酸用 p×100，扩展字母 X/Z/J/B/U/O 用该位置的背景加权均值
        for pos, row in enumerate(smatrix):
            rowmean = int(round(sum(freqs[a] * row[i] for i, a in enumerate(AA20))))
            vals = {AA20[i]: row[i] for i in range(20)}
            for L in LETTERS:
                v = vals.get(L, rowmean)
                if L == 'X':
                    v = rowmean
                mot_out.append('motif_%d@%d@%s %d' % (mid, pos, L, v))
        # 可选收紧：把 CDF 右移 delta，使阈值处的期望窗口数 ≤ target_exp
        nwin_all = sum(max(0, len(sq) - width + 1) for sq in seqs)
        delta = 0
        if a.target_exp is not None and t5 < len(tail):
            while t5 + delta < len(tail) and tail[t5 + delta] * nwin_all > a.target_exp:
                delta += 1
            if delta:
                print('  motif_%d：收紧标定 +%d 分（期望命中 %.0f → %.0f 个窗口）'
                      % (mid, delta, tail[t5] * nwin_all, tail[t5 + delta] * nwin_all), file=sys.stderr)
                # 方向很关键：**收紧 = 同样的分数对应更大的 p**，所以要把 CDF 左移
                # （cdf_new[s] = cdf_old[s - delta]），等价于"分数必须高 delta 才能达到同样的显著性"。
                # 反过来右移会放松判据——实测右移把全基因组位点从 690 推到 15754，是错误方向。
                tail = [tail[max(0, s - delta)] for s in range(len(tail))]
                t4 = min(t4 + delta, len(tail) - 1)
                t5 = min(t5 + delta, len(tail) - 1)
        for s, v in enumerate(tail):
            store_out.append('motif_%d@%d %.12g' % (mid, s, v))
        # 报告：单个随机窗口达到阈值的概率 + 全蛋白组的期望窗口数
        nwin = sum(max(0, len(s) - width + 1) for s in seqs)
        report.append(dict(id=mid, name=name, width=width,
                           thr4=t4, thr5=t5, p_at_thr5=tail[min(t5, len(tail) - 1)],
                           exp_hits_1e5=tail[t5] * nwin if t5 < len(tail) else 0.0,
                           exp_hits_1e4=tail[t4] * nwin if t4 < len(tail) else 0.0,
                           maxscore=len(tail) - 1, cdf_len=len(tail)))
        print('  motif_%d (%s)：宽 %d，p<1e-4 需得分 ≥%d，p<1e-5 需 ≥%d；'
              '全蛋白组期望命中 %.1f 个窗口（p<1e-5）'
              % (mid, name, width, t4, t5, report[-1]['exp_hits_1e5']), file=sys.stderr)

    with open(a.out_prefix + '.mot.txt', 'w') as o:
        o.write('\n'.join(mot_out) + '\n')
    with open(a.out_prefix + '.store.txt', 'w') as o:
        o.write('\n'.join(store_out) + '\n')
    with open(a.out_prefix + '.report.tsv', 'w') as o:
        o.write('motif_id\tname\twidth\tscore_p1e-4\tscore_p1e-5\texp_windows_1e5\texp_windows_1e4\tmax_score\n')
        for r in report:
            o.write('%d\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%d\n' % (
                r['id'], r['name'], r['width'], r['thr4'], r['thr5'],
                r['exp_hits_1e5'], r['exp_hits_1e4'], r['maxscore']))
    print('写出 %s.mot.txt / .store.txt / .report.tsv' % a.out_prefix)


if __name__ == '__main__':
    main()
