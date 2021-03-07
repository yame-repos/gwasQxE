# gwasQxE
This repository contains phenotype data, genotype data and an R package ("gwasQxE") corresponding to a manuscript titled "Exploring Efficient Linear Mixed Models to Detect Quantitative Trait Locus-by-environment Interactions".
The main objective of the study was to perform genome-wide association study to detect Quantitative Trait Locus-by-environment (QxE) Interactions.

The following files are included with this distribution.

  - R package:	gwasQxE
  - Genotype data:	SL96_geno_QxE.csv
  - Phenotype data:	SL96_pheno_QxE.csv
  - R source code for results in the manuscript:	simulation_note.R

<!-- end list -->


## Example

Following libraries are required to perform the analyses in note.R.
``` r
require(gaston)
require(qqman)
require(gwasQxE)
```

The simplest code to run gwasQxE is as follows
``` r 
geno = read.csv("SL96_geno_QxE.csv")
pheno = read.csv("SL96_pheno_QxE.csv")

res = gwasQxE(geno,
        pheno,
        trait = "Fruit_weight",
        Env = "season",
        n.PC = 2)
```

<!-- end list -->
