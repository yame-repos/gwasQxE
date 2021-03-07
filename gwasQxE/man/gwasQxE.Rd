\name{gwasQxE}
\alias{gwasQxE}
\title{
Genome-wide association study to detect QTL-by-Environment interaction (QxE) effect
}
\description{
This function detects QxE effect by comparing linear mixed models with assumption that a QTL has same effect in different environments and the QTL has QxE effect. This function test only additive  effect. For dominant effect, see 'gwasQxEd'.
}
\usage{
gwasQxE(geno,
        pheno,
        trait,
        Env,
        n.PC,
        Fixed,
        n.core,
        scale.pheno.by.trial)
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
\value{
Returns a data frame where the columns are the marker name, chromosome, position, -log10(p) for additive main effect (P.ame), all QTL effects (P.all), and QxE terms (P.int), additive effect (Add), additive effect in each environment (effect@~~), Wald test score for the additive effects (wald@~~), and minor allele frequency (MAF).
}
\author{
Eiji Yamamoto
}
\seealso{
gwasQxEd
}
