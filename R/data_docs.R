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

#' [DEPRECATED]
#'
#' This piece of data is deprecated and will be dismissed in next
#' release. Please use \link{GROAN.KI} instead.
#'
#' @format A data frame with 103 rows and 103 variables. Row and
#' column names are pea KI lines.
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.pea.kinship"

#' [DEPRECATED]
#'
#' This piece of data is deprecated and will be dismissed in next
#' release. Please use \link{GROAN.KI} instead.
#'
#' @format A data frame with 103 rows and 647 variables. Each
#' row represent a pea KI line, each column a SNP marker
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.pea.SNPs"

#' [DEPRECATED]
#'
#' This piece of data is deprecated and will be dismissed in next
#' release. Please use \link{GROAN.KI} instead.
#'
#' @format A named array with 103 slots.
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.pea.yield"

#' Example data for pea KI lines
#'
#' This list contains all data required to run GROAN examples. It refers to a pea experiment
#' with 103 lines coming from a biparental Kaspa x Isard cross.
#'
#' @format A list with the following fields:
#' \itemize{
#'  \item{\emph{"GROAN.KI$yield"}}: named array with 103 slots, containing data on grain yield [t/ha]
#'  \item{\emph{"GROAN.KI$SNPs"}}: data frame with 103 rows and 647 variables. Each row is
#'  a pea KI line, each column a SNP marker. Values can either be 0, 1, or 2,
#' representing the three possible genotypes (AA, Aa, and aa, respectively).
#'  \item{\emph{"GROAN.KI$kinship"}}: square dataframe containing the realized kinships between
#' all pairs of each of the 103 pea KI lines. Values were computed following the
#' \href{https://www.jstor.org/stable/25681325}{Astle & Balding metric}.
#' Higher values represent a higher degree of genetic similarity between
#' lines. This metric mainly accounts for additive genetic contributions
#' (as an alternative to dominant contributions).
#' }
#'
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.KI"

#' Example data for pea AI lines
#'
#' This list contains all data required to run GROAN examples. It refers to a pea experiment
#' with 105 lines coming from a biparental Attika x Isard cross.
#'
#' @format A list with the following fields:
#' \itemize{
#'  \item{\emph{"GROAN.AI$yield"}}: named array with 105 slots, containing data on grain yield [t/ha]
#'  \item{\emph{"GROAN.AI$SNPs"}}: data frame with 105 rows and 647 variables. Each row is
#'  a pea AI line, each column a SNP marker. Values can either be 0, 1, or 2,
#' representing the three possible genotypes (AA, Aa, and aa, respectively).
#'  \item{\emph{"GROAN.AI$kinship"}}: square dataframe containing the realized kinships between
#' all pairs of each of the 105 pea AI lines. Values were computed following the
#' \href{https://www.jstor.org/stable/25681325}{Astle & Balding metric}.
#' Higher values represent a higher degree of genetic similarity between
#' lines. This metric mainly accounts for additive genetic contributions
#' (as an alternative to dominant contributions).
#' }
#'
#' @source Annicchiarico et al., \emph{GBS-Based Genomic Selection
#' for Pea Grain Yield under Severe Terminal Drought}, The Plant Genome,
#' Volume 10. \doi{doi:10.3835/plantgenome2016.07.0072}
"GROAN.AI"
