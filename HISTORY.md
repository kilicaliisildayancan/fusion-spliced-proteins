
While looking for [proteins with multiple GPCR domains](https://github.com/kilicaliisildayancan/multi-GPCR), we had an issue with data quality. With low quality protein sequences from large protein databases, it was hard to prove that these sequences actually:
- existed at the genomic sequence level
- were transcribed at the transcriptome level 
- were produced into proteins, transported to the membrane
- were functional (if all the above is true, then I guess this is a trivial case)

The main problem was that all we found were in low-quality data sources, and the nearest high-quality source lacked them. From an evolutionary standpoint, if Wild Yak has a multi-domain GPCR protein, the close and high-quality data source Cow (the milka cows) should have it too. This happened over and over again, with no avail.

Then, searching through better-quality data, I stumbled across the relationship between two chemokine receptors that were close in proximity, CCR2 and CCR5. 
Conventionally, they are different proteins, and are translated separately for different functional purposes (have different ligands/substrates).
The most interesting thing was that, their *order* (CCR2 always upstream of CCR5), their *relative orientation* and *intergenic distance* seemed to be highly constrained across evolution. If an animal had these two receptors, they were together, facing the same direction. 

This gives rise to the possibility that the transcription machinery starts from one of them, and then transcribes until the last part of CCR5. Since transcription is not as simple as Gene A -> Transcript A -> Chain A, there are a lot of complicated machinery that changes the respective entity in the process. Specifically, splicing. We'll get to that in one second. 

The problem arises from the conservation and lack of explanation. If they *could* drift apart, change directions, mutate, or something, it is likely that *they would*. *If they don't*, there has to be an explanation. While random chance is still possible, this project aims to investigate if it is not by random chance, and if there is an important role that puts these two proteins together.

The assumption we make is that they are transcribed as one strand of mRNA. Otherwise other adaptor proteins that use the structure of the genome might be involved that require their intergenic distance to be conserved. The relative orientation still requires some sort of explanation here regardless.

There are alternatives to the reason of this phenomenon: 

1. It is done to conserve the relative concentration in the cells they are produced. In specific cases, they live in complete co-dependence.
2. There is a splicing alternative that produces a multi-domain GPCR protein with the two receptors.
3. There is a splicing alternative that combines parts of two proteins into a different GPCR.
4. It happens by random chance.
5. Some other explanation that is not listed here. 

In this research project, we're investigating the reason behind the evolutionary pressure that keeps these two genes *side by side*, and in the *same relative orientation*.

This will start out as testing to see if the second (2.) or third (3.) alternative has explanatory potential. The following search pipeline is dedicated to confirming or rejecting these hypotheses.

We are generalizing it to leverage evidence from evolution and other data sources. It is highly unlikely that this is a one-off thing (and if so, we can't call it fusion-spliced proteins, since it'd be a special case for the two chemokine receptors only), for two proteins to be translated separately and to have this topological structure (order, distance, and relative orientation).
