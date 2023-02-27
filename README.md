# Optimize potato breeding program I

This repository contains the main R scripts and empirical genotypic data that were used as a basis for the computer simulation study: "Optimal implementation of genomic selection in clone breeding programs - exemplified in potato: I. Effect of selection strategy, implementation stage, and selection intensity on short-term genetic gain" (Wu et al., 2023). These scripts were performed under R 3.6.1 and AlphaSimR 1.0.4.

Note that this manuscript is currently under revision, so changes are expected.

If you find any issues or errors, please do not hesitate to contact us:
* Po-Ya Wu (maintainer)
  - E-mail: Po-Ya.Wu@hhu.de
  - Institute for Quantitative Genetics and Genomics of Plants, Heinrich Heine University Duesseldorf 
  
* Delphine Van Inghelandt 
  - E-mail: inghelan@hhu.de
  - Institute for Quantitative Genetics and Genomics of Plants, Heinrich Heine University Duesseldorf 
  
## Data description:
* 100clones_haplotype.RData: 
  The haplotype information of 100 tetraploid potato clones used as parents in this study is stored in a list called "haplotype.all" with a length of 12 (= basic chromosome number of potato). Note that the haplotype phase was randomly determined.

* geneticmap_centromere.RData:
  The genetic map information in Morgan of 12 chromosomes as well as the estimated position of the centromere in each chromosome are stored in a list called "geno.map" and in a vector called "centromere".

* random_crossplan_list_1000rep.RData:
  The cross plan in each run/repetition was determined by randomly selecting 300 from all possible crosses in the half-diallel among the 100 parental clones, which was expressed by a matrix with $300 \times 2$. The cross plans among 1,000 repetitions were stored in a list called "crossplan.list". 
  
The other datasets were produced by performing the scripts subsequently. 
## Citation
* Wu, P.-Y., Stich, B., Renner, J., Muders, K., Prigge, V., and van Inghelandt, D. (2023). Optimal implementation of genomic selection in clone breeding programs - exemplified in potato: I. Effect of selection strategy, implementation stage, and selection intensity on short-term genetic gain. Under review .