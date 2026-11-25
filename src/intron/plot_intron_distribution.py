#!/usr/bin/env python3
"""Plot the intron-length distribution produced by intron_distribution.py.

Reads <input_dir>/{distribution.tsv, summary.tsv} and writes a histogram of
intron lengths (log-log, annotated with median / mean / P95) as PNG/SVG into
<output_dir>. Requires matplotlib (conda env `bio` on akkad).

Run from the repo root:
    python src/intron/plot_intron_distribution.py --config src/intron/plot_config.json
"""
from __future__ import annotations

import argparse
import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

VERSION = "1.0.0"


def read_tsv(path: str):
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        rows = [row for row in r if row]
    return header, rows


def read_summary(path: str) -> dict:
    _, rows = read_tsv(path)
    return {row[0]: row[1] for row in rows if len(row) >= 2}


def fmt_bp(x: float) -> str:
    if x >= 1_000_000:
        return f"{x / 1_000_000:g} Mb"
    if x >= 1_000:
        return f"{x / 1_000:g} kb"
    return f"{x:g} bp"


def main(argv=None):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Plot the intron-length distribution.")
    ap.add_argument("--config", default=os.path.join(script_dir, "plot_config.json"))
    ap.add_argument("--input-dir", help="override input_dir")
    ap.add_argument("--output-dir", help="override output_dir")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = json.load(fh)
    if args.input_dir:
        cfg["input_dir"] = args.input_dir
    if args.output_dir:
        cfg["output_dir"] = args.output_dir

    in_dir = cfg["input_dir"]
    out_dir = cfg["output_dir"]
    os.makedirs(out_dir, exist_ok=True)

    _, dist_rows = read_tsv(os.path.join(in_dir, "distribution.tsv"))
    bin_starts = [int(r[0]) for r in dist_rows]
    bin_ends = [int(r[1]) for r in dist_rows]
    counts = [int(r[3]) for r in dist_rows]

    summ = read_summary(os.path.join(in_dir, "summary.tsv"))

    def g(key):
        v = summ.get(key)
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    median = g("intron_length_median")
    mean = g("intron_length_mean")
    p95 = g("intron_length_p95")
    source = summ.get("source", "")
    build = summ.get("genome_build", "")

    fig = cfg.get("figure", {})
    figw = float(fig.get("width", 7.0))
    figh = float(fig.get("height", 5.0))
    dpi = int(fig.get("dpi", 150))

    fig_obj, ax = plt.subplots(figsize=(figw, figh), constrained_layout=True)

    widths = [e - s + 1 for s, e in zip(bin_starts, bin_ends)]
    ax.bar(bin_starts, counts, width=widths, align="edge",
           color="#4C78A8", edgecolor="none")

    if median is not None:
        ax.axvline(median, color="#C0392B", ls="--", lw=1.4,
                   label=f"median {fmt_bp(median)}")
    if mean is not None:
        ax.axvline(mean, color="#E67E22", ls=":", lw=1.4,
                   label=f"mean {fmt_bp(mean)}")
    if p95 is not None:
        ax.axvline(p95, color="#27AE60", ls="-.", lw=1.4,
                   label=f"P95 {fmt_bp(p95)}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1, 10, 100, 1_000, 10_000, 100_000, 1_000_000])
    ax.set_xticklabels(["1 bp", "10 bp", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb"])
    ax.set_xlabel("intron length (bp, log scale)")
    ax.set_ylabel("count (log)")
    ax.set_title(f"Human intron-length distribution\n{source}, {build}, protein_coding",
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)

    written = {}
    for fmt in cfg.get("formats", ["png", "svg"]):
        path = os.path.join(out_dir, f"intron_length_distribution.{fmt}")
        fig_obj.savefig(path, dpi=dpi, bbox_inches="tight")
        written[fmt] = path
    plt.close(fig_obj)

    print("=== intron distribution plot ===")
    parts = []
    if median is not None:
        parts.append(f"median {fmt_bp(median)}")
    if mean is not None:
        parts.append(f"mean {fmt_bp(mean)}")
    if p95 is not None:
        parts.append(f"P95 {fmt_bp(p95)}")
    print("  ".join(parts))
    for fmt, path in written.items():
        print(f"wrote: {path}")


if __name__ == "__main__":
    main()
