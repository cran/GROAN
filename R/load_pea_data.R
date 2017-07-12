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

# Script to store load and store data in the package
# # DATA LOADING ------------------------------------------------------------
# #phenotypes
# phenos = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/phenos.csv',
#                             stringsAsFactors = FALSE, row.names = 1)
#
# GROAN.pea.yield = phenos$grain_yield
# names(GROAN.pea.yield) = rownames(phenos)
#
# #genotypes
# GROAN.pea.SNPs = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/genos.csv',
#                           stringsAsFactors = FALSE, row.names = 1)
# #kinship
# GROAN.pea.kinship = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/kinship.AB.csv',
#                           stringsAsFactors = FALSE, row.names = 1)

# # DATA SHUFFLE AND ANONYMIZATION ------------------------------------------
# samples = names(GROAN.pea.yield)
# samples.shuff = sample(samples)
#
# GROAN.pea.yield = GROAN.pea.yield[samples.shuff]
# GROAN.pea.kinship = GROAN.pea.kinship[samples.shuff, samples.shuff]
# GROAN.pea.SNPs = GROAN.pea.SNPs[samples.shuff,]
#
# plot(GROAN.pea.yield)
#
# #new names, anonymizing everything
# samples.anon = paste(sep='_', 'sample', 1:length(GROAN.pea.yield))
# SNPs.anon = paste(sep='_', 'SNP', 1:ncol(GROAN.pea.SNPs))
#
# names(GROAN.pea.yield) = samples.anon
# rownames(GROAN.pea.kinship) = samples.anon
# colnames(GROAN.pea.kinship) = samples.anon
# rownames(GROAN.pea.SNPs) = samples.anon
# colnames(GROAN.pea.SNPs) = SNPs.anon

# # STORING INTO PACKAGE ----------------------------------------------------
# devtools::use_data(GROAN.pea.yield, GROAN.pea.SNPs, GROAN.pea.kinship)
#
# dim(GROAN.pea.kinship)
# dim(GROAN.pea.SNPs)
# length(GROAN.pea.yield)
# hist(GROAN.pea.yield)
