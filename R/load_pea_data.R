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

# Script to load and store data in the package

# # DATA LOADING KI ---------------------------------------------------------
# #phenotypes
# GROAN.pea.KI.yield = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/phenos.KI.csv',
#                               stringsAsFactors = FALSE, row.names = 1)
# #genotypes
# GROAN.pea.KI.SNPs = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/genos.KI.csv',
#                              stringsAsFactors = FALSE, row.names = 1)
# #kinship
# GROAN.pea.KI.kinship = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/kinship.AB.KI.csv',
#                                 stringsAsFactors = FALSE, row.names = 1)
#
# # DATA LOADING AI ---------------------------------------------------------
# #phenotypes
# GROAN.pea.AI.yield = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/phenos.AI.csv',
#                               stringsAsFactors = FALSE, row.names = 1)
# #genotypes
# GROAN.pea.AI.SNPs = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/genos.AI.csv',
#                              stringsAsFactors = FALSE, row.names = 1)
# #kinship
# GROAN.pea.AI.kinship = read.csv('~/research/GROAN.extra/GROAN.pea/data/to_packet/kinship.AB.AI.csv',
#                                 stringsAsFactors = FALSE, row.names = 1)
#
# # DATA SHUFFLE AND ANONYMIZATION KI ---------------------------------------
# samples.KI.shuff = sample(rownames(GROAN.pea.KI.yield))
# GROAN.pea.KI.yield = data.frame(
#   row.names = samples.KI.shuff,
#   grain_yield = GROAN.pea.KI.yield[samples.KI.shuff,]
# )
# GROAN.pea.KI.kinship = GROAN.pea.KI.kinship[samples.KI.shuff, samples.KI.shuff]
# GROAN.pea.KI.SNPs = GROAN.pea.KI.SNPs[samples.KI.shuff,]
#
# #new names, anonymizing everything
# samples.KI.anon = paste(sep='_', 'sample', 1:nrow(GROAN.pea.KI.yield))
# SNPs.KI.anon = paste(sep='_', 'SNP', 1:ncol(GROAN.pea.KI.SNPs))
#
# rownames(GROAN.pea.KI.yield) = samples.KI.anon
# rownames(GROAN.pea.KI.kinship) = samples.KI.anon
# colnames(GROAN.pea.KI.kinship) = samples.KI.anon
# rownames(GROAN.pea.KI.SNPs) = samples.KI.anon
# colnames(GROAN.pea.KI.SNPs) = SNPs.KI.anon
#
# # DATA SHUFFLE AND ANONYMIZATION AI ---------------------------------------
# samples.AI.shuff = sample(rownames(GROAN.pea.AI.yield))
# GROAN.pea.AI.yield = data.frame(
#   row.names = samples.AI.shuff,
#   grain_yield = GROAN.pea.AI.yield[samples.AI.shuff,]
# )
# GROAN.pea.AI.kinship = GROAN.pea.AI.kinship[samples.AI.shuff, samples.AI.shuff]
# GROAN.pea.AI.SNPs = GROAN.pea.AI.SNPs[samples.AI.shuff,]
#
# #new names, anonymizing everything
# samples.AI.anon = paste(sep='_', 'sample', 1:nrow(GROAN.pea.AI.yield))
# SNPs.AI.anon = paste(sep='_', 'SNP', 1:ncol(GROAN.pea.AI.SNPs))
#
# rownames(GROAN.pea.AI.yield) = samples.AI.anon
# rownames(GROAN.pea.AI.kinship) = samples.AI.anon
# colnames(GROAN.pea.AI.kinship) = samples.AI.anon
# rownames(GROAN.pea.AI.SNPs) = samples.AI.anon
# colnames(GROAN.pea.AI.SNPs) = SNPs.AI.anon
#
# # STORING INTO PACKAGE ----------------------------------------------------
# #for compatibility, KI appears twice, as several split variables and as single list
# GROAN.pea.yield = GROAN.pea.KI.yield$grain_yield
# names(GROAN.pea.yield) = rownames(GROAN.pea.KI.yield$grain_yield)
# GROAN.pea.SNPs = GROAN.pea.KI.SNPs
# GROAN.pea.kinship = GROAN.pea.KI.kinship
#
# GROAN.KI = list(
#   SNPs = GROAN.pea.KI.SNPs,
#   yield = GROAN.pea.KI.yield$grain_yield,
#   kinship = GROAN.pea.KI.kinship
# )
# names(GROAN.KI$yield) = rownames(GROAN.pea.KI.yield$grain_yield)
#
# GROAN.AI = list(
#   SNPs = GROAN.pea.AI.SNPs,
#   yield = GROAN.pea.AI.yield$grain_yield,
#   kinship = GROAN.pea.AI.kinship
# )
# names(GROAN.AI$yield) = rownames(GROAN.pea.AI.yield$grain_yield)
#
# # remember to move to GROAN root
# devtools::use_data(GROAN.pea.yield, GROAN.pea.SNPs, GROAN.pea.kinship, GROAN.AI, GROAN.KI, overwrite = TRUE)
#
# # TEST WITH GROAN ---------------------------------------------------------
# if(TRUE){
#   library(GROAN)
#
#   nds.KI.split = createNoisyDataset(name = 'KI split',
#                               genotypes = GROAN.pea.SNPs,
#                               phenotypes = GROAN.pea.yield)
#   nds.KI.list = createNoisyDataset(name = 'KI list',
#                                    genotypes = GROAN.KI$SNPs,
#                                    phenotypes = GROAN.KI$yield)
#   nds.AI.list = createNoisyDataset(name = 'AI list',
#                                    genotypes = GROAN.AI$SNPs,
#                                    phenotypes = GROAN.AI$yield)
#   wb = createWorkbench(reps = 30)
#
#   res = GROAN.run(nds.KI.split, wb)
#   res = rbind(res, GROAN.run(nds.KI.list, wb))
#   res = rbind(res, GROAN.run(nds.AI.list, wb))
#
#   plotResult(res)
# }
