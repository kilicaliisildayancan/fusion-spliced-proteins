# Finding fusion-spliced proteins

For how we got to the hypothesis, check [here](./HISTORY.md).

## What

Fusion-spliced proteins are a hypothetical species of proteins that exist as a group of tandem protein-coding genes that are transcribed together and spliced into a single functional protein with domains from both protein transcripts. The constraint here is that the individual proteins have to be translated from non-overlapping exon sequences, exist and function separately when required. Or else it's just regular alternative splicing.

Let's take two genes, A and B. Most likely, they are transcribed independently (i.e. the mRNA transcripts are not linked), and translated independently (i.e. get processed by the ribosome separately, have separate chains). This would categorize them as two separate genes, transcripts, and proteins. Here, the two genes have their own splicing alternatives (only one if no other variant exists).

Then, the fusion-spliced alternative, protein gene AB, would be some combination of A's exons, optionally the integenic region's exons, and B's exons, taking at least some part of A and some part of B. 

This is visualized in this figure:

![Overview of fusion-spliced protein genes](./figures/fusion-splicing-abstract.png)

Even in this hypothetical scenario, we can think of some constraints that focus our search methods:
- Genes face the same way
- The STOP codon in the first gene is skipped by alternative splicing, and the last part of the sequence is joined with a protein coding sequence (rather than 5' UTR or a non-coding region)
- There is no frameshift in the splicing alternative
- Intergenic distance is below the threshold for the transcription machinery
- (Trivial) Genes are independently produced to become functional proteins

## Why is it interesting?

There are two different aspects that make this interesting, one is an evolutionary aspect, and the other one is novelty. 

Assuming they exist and we find them, this demonstrates another mechanism which shows how proteins evolve over time. With eukaryotic transcription and translation mechanisms, exons can evolve independently, and be shuffled around without the requirement of independent random evolution of proteins through indels and point mutations. Two genes mix and match exons, and through time, become a single gene that has all the splicing alternatives of the union of all. 

Novelty stems from the discovery of new proteins from the existing genomes, allowing researchers to investigate more proteins, potentially uncovering key proteins in different processes.

# Methods

The premise of the fusion-splicing alternative is based on the splicing mechanism in eukaryotes. The main source of evidence will come from RNA-seq data and the analysis of that. Supporting evidence can be made with evolutionary conservation and protein structure analysis after predictions.

## Finding candidates

Although the search started with two GPCR proteins, the search pipeline aims to cast a wide net to catch all potential cases. The project will analyze RNA-seq data and will try to find fitting candidates. Then, for these cases, we plan to gather more supporting evidence from evolutionary conservation and maybe other aspects.

### Creating the list of candidates

To start looking at the genomes of species and filtering out the gene pairs that could fit our criterias, we need a threshold distance that the transcription mechanism is noted to tolerate. If two genes are spliced together, the intermediate section will be an intron. 
If we know the distribution of intron lengths, we can start finding gene pairs that are close enough to be linked together by an intron.

We take the species' genes that are alternatively spliced, and get the intron positions of each gene. From this, we will get the basepair difference between exons, to have a rough measure of how long an intron can be.

### Looking at expression profiles


### Looking for an RNA transcript that joins the CDS from one gene to the next


### Results
