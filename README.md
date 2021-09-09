
# EPIMUTESTR

EPIstasis MUTations ESTImatoR (EPIMUTESTR) is a machine learning pipeline based on Relief nearest-neighbor feature selection algorithm. It uses the advantage of a multigene interaction (epistasis) approach to find driver genes which link to different tumor types. 

#### Websites

[Lichtarge Lab](http://lichtargelab.org)

[Evolutionary Action](http://eaction.lichtargelab.org)

[Lichtarge Lab Github](https://github.com/LichtargeLab)

##### Related References

[Katsonis P, Lichtarge O., A formal perturbation equation between genotype and phenotype determines the evolutionary action of protein coding variations on fitness, Genome Research, 2014 Sep 12. pii: gr.176214.114.](https://pubmed.ncbi.nlm.nih.gov/25217195/)

[Parvandeh et al., Consensus features nested cross-validation, Bioinformatics, 2020](https://doi.org/10.1093/bioinformatics/btaa046)

### To run
1. Download the public [MC3 MAF](https://gdc.cancer.gov/about-data/publications/mc3-2017)

2. Annotate with [Evolutionary Action](http://eaction.lichtargelab.org)

3. Preprocess the EA annotated MC3 data set as input

      > Rscript 1_Preprocess_TCGA_MC3.R
      
4. Construct a case matrix (M) of samples (i) by genes (j), where each entity (Mij) is the maximum EA score among all mutations per gene. And, construct an empty control matrix (Nij) of the same size of case matrix and fill with the synthetic EA scores by maintaining the same mutational signature rate from the case matrix. 

      > Rscript 2_EPIMUTESTR_MC3.R
      
5. Generate significant genes lists

      > Rscript 3_Generate_sig_lists.R
      
6. Evaluate and compare with state-of-the-art methods 

      > Rscript 4_Criteria_for_successs
      
7. Generate the figures

      > Rscript 5_Visulatization.R
      
### Dependencies

```
install.packages(c('CORElearn', 'Rcpp', 'dplyr', 'parallel', 'foreach', 'doParallel'))

```

Other R packages

```
install.packages(c('clusterProfiler', 'org.Hs.eg.db', 'DOSE', 'dndscv', 'ggplot2', 'ggpubr'))

```

Python packages

```
from Bio import Entrez

```

