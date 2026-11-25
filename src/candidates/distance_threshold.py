#!/usr/bin/env python3
"""Find protein_coding gene pairs whose intergenic CDS gap is within a threshold.

For every protein_coding gene, take its 3'-most CDS position (the last coding
base, i.e. the 3' end in transcription direction) and pair it with every
downstream, SAME-STRAND gene whose 5'-most CDS start lies within
`intergenic_distance_threshold` bases. Overlaps (negative gaps) are kept.

Genes tagged "readthrough_gene" (already-annotated conjoined/fusion loci) are
excluded by default (config `exclude_readthrough_genes`), so the candidate list
contains only novel pairs; the excluded readthrough genes are written separately
to known_readthroughs.tsv.

Outputs (written to <output_dir>):
    candidates.tsv          every qualifying pair with CDS coordinates and gap
    known_readthroughs.tsv  protein_coding genes tagged readthrough_gene
    summary.tsv             gene/pair counts and gap statistics
    provenance.json         source, version, params, GTF hash/size

Run from the repo root:
    python src/candidates/distance_threshold.py --config src/candidates/config.json
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.request
from collections import defaultdict

VERSION = "1.1.0"


def parse_attributes(s: str) -> dict:
    """Parse a GTF attribute string ('key "value"; ...') into a dict."""
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


def parse_genes(gtf_path: str, transcript_type: str, exclude_readthrough: bool = True):
    """Return (genes, n_protein_coding_total, readthrough_genes).

    genes: gene_id -> {gene_id, gene_name, chrom, strand, cds_5start, cds_3end}
           for every protein_coding gene that has at least one CDS segment and
           (when exclude_readthrough) is not tagged "readthrough_gene".
    n_protein_coding_total: count of protein_coding gene lines (incl. no-CDS).
    readthrough_genes: gene_id -> gene_name for protein_coding genes tagged
           "readthrough_gene" (already-annotated fusion loci).
    """
    cds_segs = defaultdict(list)
    meta = {}
    coding_gene_ids = set()
    readthrough_genes = {}

    with open_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]

            if feature == "gene":
                attrs = parse_attributes(fields[8])
                if attrs.get("gene_type") == "protein_coding":
                    gid = attrs.get("gene_id")
                    coding_gene_ids.add(gid)
                    if exclude_readthrough and "readthrough_gene" in fields[8]:
                        readthrough_genes[gid] = attrs.get("gene_name", "")
            elif feature == "CDS":
                attrs = parse_attributes(fields[8])
                if attrs.get("transcript_type") != transcript_type:
                    continue
                gid = attrs.get("gene_id")
                if not gid:
                    continue
                cds_segs[gid].append((int(fields[3]), int(fields[4])))
                if gid not in meta:
                    meta[gid] = {
                        "gene_name": attrs.get("gene_name", ""),
                        "chrom": fields[0],
                        "strand": fields[6],
                    }

    genes = {}
    for gid, segs in cds_segs.items():
        if exclude_readthrough and gid in readthrough_genes:
            continue
        cds_min = min(s for s, _ in segs)
        cds_max = max(e for _, e in segs)
        strand = meta[gid]["strand"]
        if strand == "+":
            cds_5start, cds_3end = cds_min, cds_max
        else:
            cds_5start, cds_3end = cds_max, cds_min
        genes[gid] = {
            "gene_id": gid,
            "gene_name": meta[gid]["gene_name"],
            "chrom": meta[gid]["chrom"],
            "strand": strand,
            "cds_5start": cds_5start,
            "cds_3end": cds_3end,
        }

    return genes, len(coding_gene_ids), readthrough_genes


def find_pairs(genes: dict, threshold: int):
    """Yield (a, b, gap, n_between) for same-strand downstream pairs within threshold.

    gap = number of bases between a's 3' CDS end and b's 5' CDS start
          (0 = adjacent, negative = overlapping).
    """
    groups = defaultdict(list)
    for g in genes.values():
        groups[(g["chrom"], g["strand"])].append(g)

    for (chrom, strand), group in groups.items():
        # transcription order: 5' -> 3'
        if strand == "+":
            group.sort(key=lambda g: g["cds_5start"])
        else:
            group.sort(key=lambda g: g["cds_5start"], reverse=True)

        for i in range(len(group)):
            a = group[i]
            for j in range(i + 1, len(group)):
                b = group[j]
                if strand == "+":
                    gap = b["cds_5start"] - a["cds_3end"] - 1
                else:
                    gap = a["cds_3end"] - b["cds_5start"] - 1
                if gap > threshold:
                    break
                yield a, b, gap, j - i - 1


def percentile(sorted_vals, p: float):
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    k = (len(sorted_vals) - 1) * p
    lo = int(k // 1)
    hi = int(-(-k // 1))
    if lo == hi:
        return float(sorted_vals[lo])
    return sorted_vals[lo] * (hi - k) + sorted_vals[hi] * (k - lo)


def write_outputs(cfg, genes, pairs, gtf_path, readthrough_genes):
    out_dir = cfg["output_dir"]
    os.makedirs(out_dir, exist_ok=True)

    cand_path = os.path.join(out_dir, "candidates.tsv")
    rt_path = os.path.join(out_dir, "known_readthroughs.tsv")
    summary_path = os.path.join(out_dir, "summary.tsv")
    prov_path = os.path.join(out_dir, "provenance.json")

    gaps = []
    n_plus = 0
    n_minus = 0
    n_adjacent = 0
    n_skip = 0

    cols = ["gene_a_id", "gene_a_name", "gene_b_id", "gene_b_name", "chrom", "strand",
            "a_cds_5start", "a_cds_3end", "b_cds_5start", "b_cds_3end",
            "intergenic_distance", "n_between"]

    with open(cand_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for a, b, gap, n_between in pairs:
            gaps.append(gap)
            if a["strand"] == "+":
                n_plus += 1
            else:
                n_minus += 1
            if n_between == 0:
                n_adjacent += 1
            else:
                n_skip += 1
            w.writerow([a["gene_id"], a["gene_name"], b["gene_id"], b["gene_name"],
                        a["chrom"], a["strand"],
                        a["cds_5start"], a["cds_3end"], b["cds_5start"], b["cds_3end"],
                        gap, n_between])

    with open(rt_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "gene_name"])
        for gid in sorted(readthrough_genes):
            w.writerow([gid, readthrough_genes[gid]])

    gaps.sort()
    n_pairs = len(gaps)

    summary_rows = [
        ("source", cfg.get("source", "")),
        ("genome_build", cfg.get("genome_build", "")),
        ("transcript_type", cfg["transcript_type"]),
        ("intergenic_distance_threshold", cfg["intergenic_distance_threshold"]),
        ("exclude_readthrough_genes", cfg.get("exclude_readthrough_genes", True)),
        ("n_readthrough_genes_excluded", len(readthrough_genes)),
        ("n_genes_with_cds", len(genes)),
        ("n_pairs", n_pairs),
        ("n_pairs_plus_strand", n_plus),
        ("n_pairs_minus_strand", n_minus),
        ("n_adjacent_pairs", n_adjacent),
        ("n_skip_pairs", n_skip),
        ("pair_gap_min", gaps[0] if gaps else ""),
        ("pair_gap_median", percentile(gaps, 0.50) if gaps else ""),
        ("pair_gap_max", gaps[-1] if gaps else ""),
    ]
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        for k, v in summary_rows:
            w.writerow([k, v])

    sha = hashlib.sha256()
    size = os.path.getsize(gtf_path)
    with open(gtf_path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            sha.update(chunk)
    provenance = {
        "script": "src/candidates/distance_threshold.py",
        "version": VERSION,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "source": cfg.get("source", ""),
        "genome_build": cfg.get("genome_build", ""),
        "transcript_type": cfg["transcript_type"],
        "intergenic_distance_threshold": cfg["intergenic_distance_threshold"],
        "exclude_readthrough_genes": cfg.get("exclude_readthrough_genes", True),
        "gap_definition": "bases between a's 3' CDS end and b's 5' CDS start "
                          "(0 = adjacent, negative = overlap); same-strand only",
        "gtf_path": gtf_path,
        "gtf_url": cfg.get("gtf_url", ""),
        "gtf_size_bytes": size,
        "gtf_sha256": sha.hexdigest(),
    }
    with open(prov_path, "w") as fh:
        json.dump(provenance, fh, indent=2)
        fh.write("\n")

    return {
        "n_genes": len(genes),
        "n_pairs": n_pairs,
        "n_plus": n_plus,
        "n_minus": n_minus,
        "n_adjacent": n_adjacent,
        "n_skip": n_skip,
        "gap_min": gaps[0] if gaps else None,
        "gap_median": percentile(gaps, 0.50) if gaps else None,
        "gap_max": gaps[-1] if gaps else None,
        "n_readthrough": len(readthrough_genes),
    }


def main(argv=None):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Find protein_coding gene pairs within an intergenic threshold.")
    ap.add_argument("--config", default=os.path.join(script_dir, "config.json"))
    ap.add_argument("--gtf", help="override gtf_path")
    ap.add_argument("--output-dir", help="override output_dir")
    ap.add_argument("--threshold", type=int, help="override intergenic_distance_threshold")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = json.load(fh)
    if args.gtf:
        cfg["gtf_path"] = args.gtf
    if args.output_dir:
        cfg["output_dir"] = args.output_dir
    if args.threshold is not None:
        cfg["intergenic_distance_threshold"] = args.threshold

    gtf_path = ensure_gtf(cfg["gtf_path"], cfg.get("gtf_url", ""))
    print(f"[distance_threshold] GTF: {gtf_path}", file=sys.stderr)
    t0 = time.time()

    exclude_rt = cfg.get("exclude_readthrough_genes", True)
    genes, n_coding_total, readthrough_genes = parse_genes(
        gtf_path, cfg["transcript_type"], exclude_rt
    )
    pairs = list(find_pairs(genes, cfg["intergenic_distance_threshold"]))
    stats = write_outputs(cfg, genes, pairs, gtf_path, readthrough_genes)

    print("=== gene pair candidates ===")
    print(f"protein_coding genes (total): {n_coding_total}")
    print(f"readthrough genes excluded: {stats['n_readthrough']}")
    print(f"genes with CDS (screened): {stats['n_genes']}")
    print(f"candidate pairs: {stats['n_pairs']} "
          f"(+ {stats['n_plus']}, - {stats['n_minus']}; "
          f"adjacent {stats['n_adjacent']}, skip {stats['n_skip']})")
    print(f"pair gap min/median/max: {stats['gap_min']} / {stats['gap_median']} / {stats['gap_max']}")
    print(f"outputs: {cfg['output_dir']}")
    print(f"elapsed: {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
