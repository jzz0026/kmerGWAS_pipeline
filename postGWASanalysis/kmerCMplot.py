#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert of k-mer GWAS plotting from R/CMplot to interactive Plotly.
Generates an interactive Manhattan plot and QQ plot with hover showing position.

Usage:
  python3 kmerCMplot.py

Input files (in current working directory):
  - manhattan_input.txt    (five columns: Chromosome Position P kmer_id sequence)
  - threshold_5per         (optional; contains a single numeric -log10 threshold)

Outputs:
  - FT10_kmer_Manhattan.html
  - FT10_kmer_QQ.html
  - (optionally) PNGs if kaleido is available

"""
import os
import sys
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

# --- Config ---
INPUT_FILE = "manhattan_input.txt"
THRESH_FILE = "threshold_5per"
OUT_MAN_HTML = "FT10_kmer_Manhattan.html"
OUT_QQ_HTML = "FT10_kmer_QQ.html"
OUT_MAN_PNG = "FT10_kmer_Manhattan.png"
OUT_QQ_PNG = "FT10_kmer_QQ.png"
VALID_CHRS = {"1", "2", "3", "4", "5", "Chr1", "Chr2", "Chr3", "Chr4", "Chr5"}

# --- Read input ---
if not os.path.exists(INPUT_FILE):
    sys.exit(f"错误：找不到 {INPUT_FILE}，请把文件放在当前目录。")

df = pd.read_csv(INPUT_FILE, sep="\s+", header=None, names=["Chromosome", "Position", "P", "kmer_id", "sequence"], dtype=str)

# --- Clean data ---
df = df[df["Chromosome"] != "*"]
df = df[df["Chromosome"].isin(VALID_CHRS)]

# 3.3 统一去除 "Chr" 前缀并转换为数值
df["Chromosome"] = df["Chromosome"].str.replace("Chr", "")
df["Chromosome"] = pd.to_numeric(df["Chromosome"], errors="coerce")
df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
df["P"] = pd.to_numeric(df["P"], errors="coerce")
# kmer_id 和 sequence 保留为字符串
df["kmer_id"] = df["kmer_id"].astype(str)
df["sequence"] = df["sequence"].astype(str)

# 3.4 排序 (CMplot 绘图必须按顺序排列)
df = df.sort_values(["Chromosome", "Position"]).reset_index(drop=True)

# 3.5 使用 kmer_id 作为 SNP ID
df["SNP"] = df["kmer_id"]

# --- 4. 设置显著性阈值 (方案B：自动读取) ---
if os.path.exists(THRESH_FILE):
    try:
        raw_thresh = open(THRESH_FILE).read().splitlines()[0]
        thresh_log10 = float(raw_thresh.strip())
        my_threshold = 10 ** (-thresh_log10)
        print(f"Successfully loaded threshold: {thresh_log10} (P = {my_threshold})")
    except Exception:
        my_threshold = 1e-7
        print("Warning: failed to parse threshold_5per; using default 1e-7")
else:
    my_threshold = 1e-7
    print("Warning: 'threshold_5per' not found. Using default 1e-7.")

# --- 5. 绘制曼哈顿图 ---
df["neglog10P"] = -np.log10(df["P"].clip(lower=1e-323))
threshold_y = -math.log10(my_threshold)

# cumulative positions
chrom_sizes = df.groupby("Chromosome")["Position"].max().sort_index()
chrom_order = [int(x) for x in chrom_sizes.index]
offset = {}
cum = 0
for ch in chrom_order:
    offset[ch] = cum
    cum += int(chrom_sizes.loc[ch])

df["cum_pos"] = df.apply(lambda r: r["Position"] + offset[int(r["Chromosome"])], axis=1)

# Calculate tick positions at the middle of each chromosome
ticks = []
ticktext = []
for ch in chrom_order:
    ch_df = df[df["Chromosome"] == ch]
    if len(ch_df) > 0:
        mid = ch_df["cum_pos"].min() + (ch_df["cum_pos"].max() - ch_df["cum_pos"].min()) / 2
        ticks.append(mid)
        ticktext.append(str(ch))

colors = px.colors.qualitative.Dark24
color_map = {ch: colors[i % len(colors)] for i, ch in enumerate(chrom_order)}
marker_colors = df["Chromosome"].map(lambda x: color_map[int(x)])

fig = go.Figure()
fig.add_trace(go.Scattergl(
    x=df["cum_pos"],
    y=df["neglog10P"],
    mode="markers",
    marker=dict(color=marker_colors, size=5),
    text=df["SNP"],
    hovertemplate=("SNP: %{text}<br>Sequence: %{customdata[3]}<br>Chromosome: %{customdata[0]}<br>Position: %{customdata[1]}<br>p-value: %{customdata[2]:.3e}<br>-log10(p): %{y:.2f}<extra></extra>"),
    customdata=np.stack([df["Chromosome"].astype(int), df["Position"].astype(int), df["P"].astype(float), df["sequence"].astype(str)], axis=-1)
))

fig.add_shape(type="line", x0=0, x1=cum, y0=threshold_y, y1=threshold_y, line=dict(color="red", dash="dash"))

fig.update_layout(
    title="k-mer GWAS: FT10 (Threshold via Permutation Test)",
    xaxis=dict(title="Chromosome", tickmode="array", tickvals=ticks, ticktext=ticktext),
    yaxis=dict(title="-log10(p)"),
    showlegend=False,
    width=1200,
    height=600
)

fig.write_html(OUT_MAN_HTML)
print(f"Saved interactive Manhattan plot to {OUT_MAN_HTML}")
try:
    fig.write_image(OUT_MAN_PNG, scale=2)
    print(f"Also saved PNG {OUT_MAN_PNG}")
except Exception:
    print("kaleido not available or failed; skipping PNG export for Manhattan.")

# --- 6. 绘制 QQ 图 ---
sorted_p = np.sort(df["P"].values)
n = len(sorted_p)
obs = -np.log10(sorted_p.clip(min=1e-323))
exp = -np.log10((np.arange(1, n + 1) / (n + 1.0)))

# Keep mapping to SNPs and sequences for hover: sorted indices
sorted_idx = np.argsort(df["P"].values)
qq_snp = df.iloc[sorted_idx]["SNP"].values
qq_seq = df.iloc[sorted_idx]["sequence"].values

qq_fig = go.Figure()
qq_fig.add_trace(go.Scattergl(
    x=exp,
    y=obs,
    mode="markers",
    marker=dict(size=5, color="blue"),
    text=qq_snp,
    customdata=qq_seq,
    hovertemplate=("SNP: %{text}<br>Sequence: %{customdata}<br>Expected -log10(p): %{x:.2f}<br>Observed -log10(p): %{y:.2f}<extra></extra>")
))
maxv = max(exp.max(), obs.max())
qq_fig.add_trace(go.Scatter(x=[0, maxv], y=[0, maxv], mode="lines", line=dict(color="gray", dash="dash"), showlegend=False))
qq_fig.update_layout(title="QQ Plot of FT10 k-mer GWAS", xaxis_title="Expected -log10(p)", yaxis_title="Observed -log10(p)", width=800, height=600)

qq_fig.write_html(OUT_QQ_HTML)
print(f"Saved interactive QQ plot to {OUT_QQ_HTML}")
try:
    qq_fig.write_image(OUT_QQ_PNG, scale=2)
    print(f"Also saved PNG {OUT_QQ_PNG}")
except Exception:
    print("kaleido not available or failed; skipping PNG export for QQ.")

print("Done.")
