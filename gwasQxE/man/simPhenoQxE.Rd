\name{simPhenoQxE}
\alias{simPhenoQxE}
\title{
Generator of simulated phenotypic values
}
\description{
This function generates simulated phenotypic values using genotype data in the argument 'geno'.
}
\usage{
simPhenoQxE(geno,
            n.qtl,
            seed,
            q.Effs,
            Env,
            w_ratio,
            PVE.qtl,
            h2,
            domEff = FALSE)
}
\arguments{
\item{geno}{
Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position (either bp or cM), respectively. Columns four and higher contain the marker scores for each line, coded as (-1,0,1) = (aa,Aa,AA). Missing (NA) values are not allowed.
}
\item{n.qtl}{
Number of major QTLs in the simulated phenotype.
}
\item{seed}{
Random number generator state.
}
\item{q.Effs}{
Effect size of the major QTLs for each experimental trials. Namely, the number of elements in q.Effs must be equal to total number of experimental trials in the simulated phenotype.
}
\item{Env}{
Identifier of environmental conditions for each q.Effs value. The number of elements in Env must be equal to the number of elements in q.Effs.
}
\item{w_ratio}{
Ratio of contribution of background random effects. This must contain three values. First to third for not-specific-to-environments, genotype-by-environment, and genotype-by-trial effects, respectively.
}
\item{PVE.qtl}{
Proportion of variance explained by each major QTL in total phenotypic values (namely, not specific to environments or trials).
}
\item{h2}{
Heritability for total phenotypic values (namely, not specific to environments or trials).
}
\item{domEff}{
If TRUE, the major QTLs include dominant effects. The dominant effect is coded as (1,3,2) = (aa,Aa,AA). Default FALSE.
}
}
\value{
$pheno is data frame with the line names in the first column. The second and third columns contain  identifier for environments and experimental trials, respectively. The fourth column is the simulated phenotypic values. $simQTLs is the marker names of the simulated major QTLs. $realized_PVE.qtl is the realized PVE of each major QTL. $realized_h2 is the realized h2. $PVE.each contains realized PVE of the major QTLs in each trial.
}
\author{
Eiji Yamamoto
}
\references{
Yamamoto and Matsunaga. 2021. Exploring efficient linear mixed models to detect quantitative trait locus-by-environ- ment interactions. G3 (Bethesda) doi: 10.1093/g3journal/jkab119.
}
\seealso{
gwasQxE
}
