\name{gwasQxE_all}
\alias{gwasQxE_all}
\title{
Comparison of power to detect QTL-by-Environment interactions (QxE) between the linear mixed models.
}
\description{
This function calculate -log10(p) from the linear mixed models that consist of
various combinations of fixed genetic effect and background random genetic effect terms.
}
\usage{
gwasQxE_all(geno,
            pheno,
            trait,
            Env,
            n.PC,
            Fixed,
            n.core,
            scale.pheno.by.trial = TRUE)
}
\arguments{
\item{geno}{
Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position (either bp or cM), respectively. Columns four and higher contain the marker scores for each line, coded as (-1,0,1) = (aa,Aa,AA). Missing (NA) values are not allowed. The column names must match the line names in the "pheno" data frame.
}
\item{pheno}{
Data frame where the first column is the line name (line_name) and the second colum is information on trials. The remaining columns include environments, phenotypic values, and fixed covariate(s).
}
\item{trait}{
The name of the column that is used as phenotypic values.
}
\item{Env}{
The name of the column that is used as environmental covariate.
}
\item{n.PC}{
Number of principal components to be include as fixed effects.
}
\item{Fixed}{
The name(s) of the column(s) that should be included as fixed effect(s).
}
\item{n.core}{
Number of cores to be used for the calculation.
}
\item{scale.pheno.by.trial}{
Whether the phenotypic values are scaled for each trial or not. Default is TRUE.
}
}
\details{
This function performs GWAS using all the linear mixed models presented in the reference manuscript. The main objective of this function is a comparison of power between the linear mixed models. The best combination of the liner mixed models to detect QxE was selected based on results from this function and has been implemented in the another function "gwasQxE".
}
\value{
A data frame where the first to the third columns are the marker name, chromosome and position, respectively. The subsequent colmuns are -log10(p) values for each linear mixed model.
}
\author{
Eiji Yamamoto
}
\seealso{
gwasQxEd_all
}
