# References — readthrough transcription and chimeric / fusion RNAs

Background reading on the readthrough / cis-SAGe / chimeric-RNA phenomenon that
underlies this project (two adjacent same-strand genes transcribed as one and
spliced into a single transcript). Grouped by theme.

## Foundational (genome-wide discovery)

1. Akiva P, Toporik A, Edelheit S, Peretz Y, Diber A, Shemesh R, Novik A, Sorek R.
   Transcription-mediated gene fusion in the human genome. Genome Res. 2006;16(1):30-36.
   PMID 16344562. doi:10.1101/gr.4137606.
   Landmark: systematically identified 200+ cases of intergenic splicing (421 genes) and
   validated ~half in human tissues; proposes transcription-induced chimerism followed by
   retroposition as an evolutionary route to fused genes.

2. Parra G, Reymond A, Dabbouseh N, Dermitzakis ET, Castelo R, Thomson TM, Antonarakis SE,
   Guigo R. Tandem chimerism as a means to increase protein complexity in the human genome.
   Genome Res. 2006;16(1):37-44. PMID 16344564.
   Estimates 4-5% of tandem gene pairs can be co-transcribed into a single RNA encoding a
   putative chimeric protein — the framing closest to this project's hypothesis.

## Mechanism and regulation

3. Chanfreau GF. The Goldilocks state of transcriptional readthrough for chimeric splicing.
   Transcription. 2026;17(2):87-99. doi:10.1080/21541264.2026.2642485. PMID 41962124.
   Current review of how chimeric splicing arises and what tilts readthrough toward (or away
   from) actually producing chimeric mRNAs.

4. Grosso AR, Leite AP, Carvalho S, Matos MR, Martins FB, Vitor AC, Desterro JM,
   Carmo-Fonseca M, de Almeida SF. Pervasive transcription read-through promotes aberrant
   expression of oncogenes and RNA chimeras in renal carcinoma. eLife. 2015;4:e09214.
   doi:10.7554/eLife.09214. PMID 26575290.
   Readthrough as a source of RNA chimeras in ccRCC; junction pattern (penultimate exon of the
   5' gene to the second exon of the 3' gene) and a recurrent CTSC-RAB38 chimera.

5. Qin F, Song Z, Babiceanu M, Song Y, Facemire L, Singh R, Adli M, Li H. Discovery of
   CTCF-sensitive cis-spliced fusion RNAs between adjacent genes in human prostate cells.
   PLoS Genet. 2015;11(2):e1005001. doi:10.1371/journal.pgen.1005001.
   First genome-wide cis-SAGe screen (via CTCF knockdown); codifies the "rules" for cis-SAGe:
   same-strand neighbors <30 kb, 5' gene active, penultimate-exon to second-exon junction.

6. Modestov A, Zakharova G, Poddubskaya E, Buzdin A. Regulatory mechanisms and functional roles
   of readthrough transcripts in tumorigenesis. Int J Mol Sci. 2025;26(20):9975.
   doi:10.3390/ijms26209975. PMID 41155268.
   Broad recent review of readthrough / DoG transcripts and cis-SAGe, including the five
   structural types and the estimate that 10-20% of genes in healthy tissues produce cis-SAGe
   transcripts (up to ~55% under stress).

## Landscape in normal tissues (chimeras are not cancer-specific)

7. Nacu S, Yuan W, Kan Z, Bhatt D, Rivers CS, Stinson J, Peters BA, Modrusan Z, Jung K,
   Seshagiri S, Wu TD. Deep RNA sequencing analysis of readthrough gene fusions in human
   prostate adenocarcinoma and reference samples. BMC Med Genomics. 2011;4:11.
   doi:10.1186/1755-8794-4-11. PMID 21261984.
   RNA-seq catalog of 339 transcription-induced chimeras (TICs); expression tracks the upstream
   gene, not the downstream one.

8. Babiceanu M, Qin F, Xie Z, Jia Y, Lopez K, Janus N, Facemire L, Kumar S, et al. Recurrent
   chimeric fusion RNAs in non-cancer tissues and cells. Nucleic Acids Res. 2016;44(6):2859-2872.
   Chimeric RNAs recur in normal (non-cancer) cells — grounds them as general physiology rather
   than a cancer artifact.

9. Singh S, Qin F, Kumar S, et al. The landscape of chimeric RNAs in non-diseased tissues and
   cells. Nucleic Acids Res. 2020;48(4):1764-1778.
   Later, larger survey of chimeric RNAs across normal tissues.

## Function of chimeric proteins / RNAs

10. Frenkel-Morgenstern M, Lacroix V, Ezkurdia I, et al. Chimeras taking shape: potential
    functions of proteins encoded by chimeric RNA transcripts. Genome Res. 2012;22(7):1231-1242.
    doi:10.1101/gr.130062.111. PMID 22588898.
    Whether chimeric RNAs are translated, and what the fused proteins might do.

11. Yuan H, Qin F, Movassagh M, Park H, Golden W, Xie Z, Zhang P, Sklar J, Li H. A chimeric RNA
    characteristic of rhabdomyosarcoma in normal myogenesis process. Cancer Discov.
    2013;3(12):1394-1403. PMID 24089019.
    A chimeric RNA found in tumors is also produced during normal myogenesis — chimeras as
    normal developmental events, not just oncogenic ones.

12. Qin F, Zhang Y, Liu J, Li H. SLC45A3-ELK4 functions as a long non-coding chimeric RNA.
    Cancer Lett. 2017;404:53-61. PMID 28716526.
    A cis-SAGe fusion that encodes the same protein as one parent but acts as a lncRNA —
    chimeras need not function through a fused protein.

## Review

13. Li Z, Qin F, Li H. Chimeric RNAs and their implications in cancer. Curr Opin Genet Dev.
    2018;48:36-43. doi:10.1016/j.gde.2017.10.002. PMID 29100211.
    Terminology (cis-SAGe vs trans-splicing vs gene fusion) plus mechanisms and cancer
    implications; good entry point for the field's vocabulary.

## Notes

- These are the "known cases" that GENCODE's readthrough_gene tag corresponds to; this project
  asks the stricter question of whether the two halves also function independently and whether
  the fusion is translated and functional.
- Recurring junction pattern across independent studies: penultimate exon of the 5' gene spliced
  to the second exon of the 3' gene (see refs 4, 5) — relevant when we check candidate pairs at
  transcript level.
