#!/usr/bin/env python3
"""Genome-wide intron-length distribution from a GENCODE GTF (protein_coding).

Parses every protein_coding transcript's exons, derives introns as the genomic
gap between consecutive exons in transcription order, and writes a re-derivable
raw table plus derived lengths / summary / distribution / provenance files.

Outputs (written to <output_dir>):
    raw.tsv          every exon and intron, numbered in transcription order,
                     with genomic location and length (source of truth)
    lengths.tsv      same features, length only
    summary.tsv      counts + intron/exon length summary statistics
    distribution.tsv intron-length histogram (log2 bins)
    provenance.json  source, version, params, input-file hash/size

Run from the repo root:
    python src/intron/intron_distribution.py --config src/intron/config.json
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
import time
import urllib.request
from collections import defaultdict

VERSION = "1.0.0"


def parse_attributes(s: str) -> dict:
    """Parse a GTF attribute string ('key "value"; key2 "value2";') into a dict."""
    attrs = {}
    for token in s.split(";"):
        token = token.strip()
        if not token:
            continue
        if '"' in token:
            key, _, rest = token.partition('"')
            val = rest.rsplit('"', 1)[0]
            attrs[key.strip()] = val
        else:
            parts = token.split(None, 1)
            if parts:
                attrs[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return attrs


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def download(url: str, dest: str) -> None:
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    print(f"[download] {url} -> {dest}", file=sys.stderr)
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)


def ensure_gtf(gtf_path: str, gtf_url: str) -> str:
    if os.path.exists(gtf_path):
        return gtf_path
    if not gtf_url:
        raise SystemExit(f"GTF not found: {gtf_path} (and no gtf_url to download)")
    download(gtf_url, gtf_path)
    return gtf_path


def order_exons(exons, strand: str):
    """Sort exons (exon_number, start, end) into transcription order.

    Prefers the GTF exon_number; falls back to genomic order (descending on the
    minus strand) when exon_number is absent.
    """
    if all(e[0] and e[0] > 0 for e in exons):
        return sorted(exons, key=lambda e: e[0])
    return sorted(exons, key=lambda e: e[1], reverse=(strand == "-"))


def parse_gtf(gtf_path: str, transcript_type: str):
    """Return (meta, exons_by_tid, introns_by_tid) for the given transcript_type.

    meta:            tid -> dict(gene_id, gene_name, chrom, strand)
    exons_by_tid:    tid -> [(exon_number, start, end), ...] transcription order
    introns_by_tid:  tid -> [(start, end), ...] transcription order
    """
    meta = {}
    coding_tids = set()
    exons_by_tid = defaultdict(list)

    with open_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = parse_attributes(fields[8])

            if feature == "transcript":
                if attrs.get("transcript_type") == transcript_type:
                    tid = attrs.get("transcript_id")
                    if tid:
                        coding_tids.add(tid)
                        meta[tid] = {
                            "gene_id": attrs.get("gene_id", ""),
                            "gene_name": attrs.get("gene_name", ""),
                            "chrom": fields[0],
                            "strand": fields[6],
                        }
            elif feature == "exon":
                tid = attrs.get("transcript_id")
                if not tid:
                    continue
                if attrs.get("transcript_type") == transcript_type or tid in coding_tids:
                    exon_number = attrs.get("exon_number")
                    exon_number = (
                        int(exon_number) if (exon_number and exon_number.isdigit()) else 0
                    )
                    exons_by_tid[tid].append((exon_number, int(fields[3]), int(fields[4])))
                    if tid not in meta:
                        meta[tid] = {
                            "gene_id": attrs.get("gene_id", ""),
                            "gene_name": attrs.get("gene_name", ""),
                            "chrom": fields[0],
                            "strand": fields[6],
                        }

    introns_by_tid = {}
    for tid, exons in exons_by_tid.items():
        strand = meta[tid]["strand"]
        ordered = order_exons(exons, strand)
        exons_by_tid[tid] = ordered
        introns = []
        for i in range(len(ordered) - 1):
            _, s_a, e_a = ordered[i]
            _, s_b, e_b = ordered[i + 1]
            # The intron lies between the 3' end of exon i and the 5' start of
            # exon i+1. Which coordinate is 3'/5' depends on the strand.
            inner_a = e_a if strand == "+" else s_a
            inner_b = s_b if strand == "+" else e_b
            istart = min(inner_a, inner_b) + 1
            iend = max(inner_a, inner_b) - 1
            introns.append((istart, iend))
        introns_by_tid[tid] = introns

    return meta, exons_by_tid, introns_by_tid


def percentile(sorted_vals, p: float):
    """Linear-interpolated percentile of a sorted list (0 <= p <= 1)."""
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    k = (len(sorted_vals) - 1) * p
    lo = int(math.floor(k))
    hi = int(math.ceil(k))
    if lo == hi:
        return float(sorted_vals[lo])
    return sorted_vals[lo] * (hi - k) + sorted_vals[hi] * (k - lo)


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return round(v, 3)
    return v


def write_outputs(cfg, meta, exons_by_tid, introns_by_tid, gtf_path):
    out_dir = cfg["output_dir"]
    os.makedirs(out_dir, exist_ok=True)

    raw_path = os.path.join(out_dir, "raw.tsv")
    lengths_path = os.path.join(out_dir, "lengths.tsv")
    summary_path = os.path.join(out_dir, "summary.tsv")
    dist_path = os.path.join(out_dir, "distribution.tsv")
    prov_path = os.path.join(out_dir, "provenance.json")

    intron_lengths = []
    exon_lengths = []
    gene_ids = set()
    n_exons = 0
    n_introns = 0
    tids = sorted(exons_by_tid.keys())

    raw_cols = ["gene_id", "gene_name", "transcript_id", "chrom", "strand",
                "feature", "index", "start", "end", "length"]
    len_cols = ["gene_id", "gene_name", "transcript_id", "chrom", "strand",
                "feature", "index", "length"]

    with open(raw_path, "w", newline="") as raw_fh, open(lengths_path, "w", newline="") as len_fh:
        raw_w = csv.writer(raw_fh, delimiter="\t", lineterminator="\n")
        len_w = csv.writer(len_fh, delimiter="\t", lineterminator="\n")
        raw_w.writerow(raw_cols)
        len_w.writerow(len_cols)

        for tid in tids:
            m = meta[tid]
            gene_id = m["gene_id"]
            gene_name = m["gene_name"]
            chrom = m["chrom"]
            strand = m["strand"]
            gene_ids.add(gene_id)
            exons = exons_by_tid[tid]
            introns = introns_by_tid[tid]

            # Emit interleaved: exon1, intron1, exon2, intron2, ..., exonN.
            for i in range(len(exons)):
                _, s, e = exons[i]
                elen = e - s + 1
                exon_lengths.append(elen)
                n_exons += 1
                raw_w.writerow([gene_id, gene_name, tid, chrom, strand, "exon", i + 1, s, e, elen])
                len_w.writerow([gene_id, gene_name, tid, chrom, strand, "exon", i + 1, elen])

                if i < len(introns):
                    is_, ie = introns[i]
                    ilen = ie - is_ + 1
                    intron_lengths.append(ilen)
                    n_introns += 1
                    raw_w.writerow([gene_id, gene_name, tid, chrom, strand, "intron", i + 1, is_, ie, ilen])
                    len_w.writerow([gene_id, gene_name, tid, chrom, strand, "intron", i + 1, ilen])

    intron_lengths.sort()
    exon_lengths.sort()

    def q(vals, p):
        return percentile(vals, p)

    summary_rows = [
        ("source", cfg.get("source", "")),
        ("genome_build", cfg.get("genome_build", "")),
        ("transcript_type", cfg["transcript_type"]),
        ("n_genes", len(gene_ids)),
        ("n_transcripts", len(tids)),
        ("n_exons", n_exons),
        ("n_introns", n_introns),
        ("intron_length_min", intron_lengths[0] if intron_lengths else ""),
        ("intron_length_p25", q(intron_lengths, 0.25)),
        ("intron_length_median", q(intron_lengths, 0.50)),
        ("intron_length_mean", sum(intron_lengths) / len(intron_lengths) if intron_lengths else ""),
        ("intron_length_p75", q(intron_lengths, 0.75)),
        ("intron_length_p90", q(intron_lengths, 0.90)),
        ("intron_length_p95", q(intron_lengths, 0.95)),
        ("intron_length_p99", q(intron_lengths, 0.99)),
        ("intron_length_max", intron_lengths[-1] if intron_lengths else ""),
        ("exon_length_min", exon_lengths[0] if exon_lengths else ""),
        ("exon_length_median", q(exon_lengths, 0.50)),
        ("exon_length_mean", sum(exon_lengths) / len(exon_lengths) if exon_lengths else ""),
        ("exon_length_max", exon_lengths[-1] if exon_lengths else ""),
    ]
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        for k, v in summary_rows:
            w.writerow([k, _fmt(v)])

    bins = defaultdict(int)
    for length in intron_lengths:
        if length <= 0:
            continue
        bins[int(math.floor(math.log2(length)))] += 1
    total = sum(bins.values())
    with open(dist_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["bin_start", "bin_end", "bin_label", "count", "fraction", "cumulative_fraction"])
        cum = 0.0
        for k in sorted(bins):
            bstart = 2 ** k
            bend = 2 ** (k + 1) - 1
            c = bins[k]
            frac = c / total if total else 0.0
            cum += frac
            w.writerow([bstart, bend, f"2^{k} ({bstart}-{bend})", c, round(frac, 6), round(cum, 6)])

    sha = hashlib.sha256()
    size = os.path.getsize(gtf_path)
    with open(gtf_path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            sha.update(chunk)
    provenance = {
        "script": "src/intron/intron_distribution.py",
        "version": VERSION,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "source": cfg.get("source", ""),
        "genome_build": cfg.get("genome_build", ""),
        "transcript_type": cfg["transcript_type"],
        "gtf_path": gtf_path,
        "gtf_url": cfg.get("gtf_url", ""),
        "gtf_size_bytes": size,
        "gtf_sha256": sha.hexdigest(),
        "parameters": {"histogram": "log2 bins", "min_intron_length": 1},
    }
    with open(prov_path, "w") as fh:
        json.dump(provenance, fh, indent=2)
        fh.write("\n")

    return {
        "n_genes": len(gene_ids),
        "n_transcripts": len(tids),
        "n_exons": n_exons,
        "n_introns": n_introns,
        "intron_median": q(intron_lengths, 0.50),
        "intron_mean": sum(intron_lengths) / len(intron_lengths) if intron_lengths else None,
        "intron_max": intron_lengths[-1] if intron_lengths else None,
    }


def main(argv=None):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Genome-wide intron-length distribution from a GENCODE GTF.")
    ap.add_argument("--config", default=os.path.join(script_dir, "config.json"),
                    help="path to config.json (default: alongside this script)")
    ap.add_argument("--gtf", help="override gtf_path from the config")
    ap.add_argument("--output-dir", help="override output_dir from the config")
    ap.add_argument("--transcript-type", help="override transcript_type from the config")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = json.load(fh)
    if args.gtf:
        cfg["gtf_path"] = args.gtf
    if args.output_dir:
        cfg["output_dir"] = args.output_dir
    if args.transcript_type:
        cfg["transcript_type"] = args.transcript_type

    gtf_path = ensure_gtf(cfg["gtf_path"], cfg.get("gtf_url", ""))
    print(f"[intron_distribution] GTF: {gtf_path}", file=sys.stderr)
    t0 = time.time()

    meta, exons_by_tid, introns_by_tid = parse_gtf(gtf_path, cfg["transcript_type"])
    stats = write_outputs(cfg, meta, exons_by_tid, introns_by_tid, gtf_path)

    print("=== intron distribution ===")
    for k, v in stats.items():
        print(f"{k}: {_fmt(v)}")
    print(f"outputs: {cfg['output_dir']}")
    print(f"elapsed: {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
