# Copyright 2017 Nelson Nazzicari
# This file is part of GROAN.
#
# GROAN is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GROAN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GROAN  If not, see <http://www.gnu.org/licenses/>.

#' Data - kinship among pea lines
#'
#' A square dataframe containing the realized kinships between
#' pairs of pea lines. Values were computed following the
#' \href{http://www.jstor.org/stable/25681325}{Astle & Balding metric}.
#' Higher values represent a higher degree of genetic similarity between
#' lines. This metric mainly accounts for additive genetic contributions
#' (as an alternative to dominant contributions).
#'
#' @format A data frame with 103 rows and 103 variables. Row and
#' column names are pea lines.
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.pea.kinship"

#' Data - pea SNP genotypes
#'
#' SNP genotypes for 103 pea lines. Each line of this data frame represent
#' a single pea line. Each column a SNP marker. Values can either be 0, 1, or 2,
#' representing the three possible genotypes (AA, Aa, and aa, respectively).
#'
#' @format A data frame with 103 rows and 647 variables. Each
#' row represent a pea line, each column a SNP marker
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.pea.SNPs"

#' Data - pea phenotypes (yield)
#'
#' Phenotype dataset, containing data on pea grain yield [t/ha]
#' for each pea line.
#'
#' @format A named array with 103 slots.
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.pea.yield"
