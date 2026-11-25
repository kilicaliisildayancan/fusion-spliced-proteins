# Gene-pair candidates

Second workstream: enumerate protein_coding gene pairs close enough to be
joined by an intron, per the fusion-splicing premise.

## Method

- Source: GENCODE v50 basic, human GRCh38, `transcript_type == protein_coding`.
- For each protein_coding gene, the CDS is collapsed from all its protein_coding
  transcripts: `cds_5start` = 5'-most coding base, `cds_3end` = 3'-most coding
  base (strand-aware).
- A gene A pairs with every downstream, SAME-STRAND gene B whose 5'-most CDS
  start lies within `intergenic_distance_threshold` bases of A's 3'-most CDS
  end. Overlaps (negative gaps) are kept. Intervening genes do not block:
  a->c qualifies even with b in between (and then a->b qualifies trivially).
- Gap = number of bases between A's 3' CDS end and B's 5' CDS start
  (0 = adjacent, negative = overlapping). Same-strand only — a fused transcript
  requires the two genes to be transcribed by the same machinery.

## Threshold

`intergenic_distance_threshold` is a config variable. Default 23603 = the P95
of the intron-length distribution from `src/intron/` (95% of human protein_coding
introns are <= 23603 bp).

## Run

    python src/candidates/gene_pairs.py --config src/candidates/config.json

If the GTF in `gtf_path` is missing it is downloaded from `gtf_url`.

## Outputs

- `candidates.tsv` — each qualifying pair. Columns: gene_a_id, gene_a_name,
  gene_b_id, gene_b_name, chrom, strand, a_cds_5start, a_cds_3end,
  b_cds_5start, b_cds_3end, intergenic_distance, n_between.
- `summary.tsv` — gene/pair counts and gap statistics (metric, value).
- `provenance.json` — source, version, params, GTF path/url/hash/size.

## Annotation note (basic vs comprehensive)

Uses the basic annotation (representative transcripts, no readthroughs). CDS
extents reflect primary isoforms only; a gene with a distal alternative coding
isoform present only in the comprehensive set would have a slightly shorter CDS
here, which can overstate a gap. Using comprehensive would capture maximal CDS
extents but would also import readthrough transcripts that already span adjacent
genes (which would show as overlaps). Basic is the conservative choice and is
consistent with the intron-length threshold; swap `gtf_path` to the
comprehensive GTF to test the sensitivity.
