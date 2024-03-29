% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_docs.R
\docType{data}
\name{GROAN.AI}
\alias{GROAN.AI}
\title{Example data for pea AI lines}
\format{
A list with the following fields:
\itemize{
 \item{\emph{"GROAN.AI$yield"}}: named array with 105 slots, containing data on grain yield [t/ha]
 \item{\emph{"GROAN.AI$SNPs"}}: data frame with 105 rows and 647 variables. Each row is
 a pea AI line, each column a SNP marker. Values can either be 0, 1, or 2,
representing the three possible genotypes (AA, Aa, and aa, respectively).
 \item{\emph{"GROAN.AI$kinship"}}: square dataframe containing the realized kinships between
all pairs of each of the 105 pea AI lines. Values were computed following the
\href{https://www.jstor.org/stable/25681325}{Astle & Balding metric}.
Higher values represent a higher degree of genetic similarity between
lines. This metric mainly accounts for additive genetic contributions
(as an alternative to dominant contributions).
}
}
\source{
Annicchiarico et al., \emph{GBS-Based Genomic Selection
for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
}
\usage{
GROAN.AI
}
\description{
This list contains all data required to run GROAN examples. It refers to a pea experiment
with 105 lines coming from a biparental Attika x Isard cross.
}
\keyword{datasets}
