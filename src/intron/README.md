# Intron-length distribution

First workstream: genome-wide human intron-length distribution, to establish how
long an intron can be — the threshold used later to filter gene-pair candidates.

## Method

- Source: GENCODE v50 basic, human GRCh38, `transcript_type == protein_coding`.
- An intron is the genomic gap between two consecutive exons of the same
  transcript, in transcription order. Strand-aware: the gap lies between the
  3' end of exon i and the 5' start of exon i+1, so
  `[min(inner_a, inner_b)+1, max(inner_a, inner_b)-1]` where the inner edges are
  `end`/`start` on the plus strand and `start`/`end` on the minus strand.
- Coordinates are 1-based inclusive (GTF convention). Intron length =
  `end - start + 1`.

## Run

    python src/intron/intron_distribution.py --config src/intron/config.json

If the GTF in `gtf_path` is missing it is downloaded from `gtf_url`.
Outputs go to `results/intron/intron_distribution/` (relative to the repo root).

## Outputs

- `raw.tsv` — every exon and intron, numbered in transcription order, with
  genomic location and length. Columns: gene_id, gene_name, transcript_id,
  chrom, strand, feature (exon|intron), index, start, end, length. Source of
  truth; all other outputs can be re-derived from it.
- `lengths.tsv` — same features, length only (gene_id, gene_name,
  transcript_id, chrom, strand, feature, index, length).
- `summary.tsv` — counts + intron/exon length summary statistics
  (metric, value).
- `distribution.tsv` — intron-length histogram, log2 bins (bin_start, bin_end,
  bin_label, count, fraction, cumulative_fraction).
- `provenance.json` — source, version, params, GTF path/url/hash/size.
