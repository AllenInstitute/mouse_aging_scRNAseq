## Description/Goal:
## Running MAST at the supertype level in non-neuronals to find age-DE genes

library(MAST)
library(scrattch.hicat)
library(dplyr)
library(ggplot2)
library(data.table)
library(bigstatsr)
library(scrattch.bigcat)


################################################################################
## Load data
################################################################################
## Cleaned metadata data file
load("anno.df.clean_cluster0916_v231221_20240206.rda")

## Loading cell by gene matrix in dgC sparse matrix format (not included in github)
load("norm.dat_objects/nonneuronal-imn/norm.dat.rda")


################################################################################
## Deciding which subclasses to include for analysis
# ################################################################################
anno.df.select = anno.df.clean[grepl("NN", anno.df.clean$supertype_label),]
tmp = as.data.frame(table(anno.df.select[,c("supertype_label", "age_cat")]))
tmp = tmp[tmp$Freq > 100,]
supertype.keep = names(table(tmp$supertype_label)[table(tmp$supertype_label) == 2])
common = intersect(colnames(norm.dat), anno.df.clean$sample_id)
anno.df.select = anno.df.clean[anno.df.clean$sample_id %in% common,]


################################################################################
## MAST per supertype (using parallel function)
################################################################################
source("functions/degenes_age_mast_simplified_qcscore.R") ## custom function for running MAST in parallel
setwd("MAST_nnsupertype/")

de.gene.results = get_age_degenes_mast_p(anno.df.select, norm.dat, groups = supertype.keep, split.by = "supertype_label", save.as.tmp = F,
                                             maxcells = 2000, freq_expressed = 0.1, mc.cores = 10)

saveRDS(de.gene.results, file = "de.gene.result_cluster0916_v231221_nnsupertype.rds")


################################################################################
## Assemble MAST results
################################################################################
opath = "MAST_nnsupertype"

## Load results
load("anno.df.clean_cluster0916_v231221_20240206.rda")
results = readRDS("MAST_nnsupertype/de.gene.result_cluster0916_v231221_nnsupertype.rds")

## Exclude genes
load("r_objects/sex.genes.rda") #sex-biased genes
load("r_objects/genotype.genes.rda") #genotype-biased genes
load("r_objects/iegenes.rda") #immediate early genes
genes.exclude = c(iegenes.rapidprg, iegenes.delayedprg, sex.genes, genotype.genes)

## Cutoffs for age-DE genes
fcthreshold = 1
pthreshold = 0.01

## Assemble results
de.gene.results = lapply(results, function(x) x[["result"]])

for(i in names(de.gene.results)){
  de.gene.results[[i]]$subclass = i
}
de.results.long0 = do.call(rbind, de.gene.results)
saveRDS(de.results.long0, file = paste0(opath, "de.results.long0_allresults.rds"))

de.gene.results = lapply(de.gene.results, function(x) x[x$padjust < pthreshold & abs(x$coef) > fcthreshold,])
de.gene.results = lapply(de.gene.results, function(x) x[order(x$coef, decreasing = T),])
de.gene.results = lapply(de.gene.results, function(x) x[!x$primerid %in% genes.exclude,]) ## excluding sex, genotype, iegenes
sort(unlist(lapply(de.gene.results, nrow)))
subclass.results = as.data.frame(sort(unlist(lapply(de.gene.results, nrow))))
names(subclass.results) = "de.gene.count"
subclass.results$subclass = row.names(subclass.results)

tmp = unlist(lapply(results, function(x) length(x$select.cells)))
names(tmp) = names(results)
subclass.results$ncell.tested = tmp[subclass.results$subclass]

saveRDS(subclass.results, file = paste0(opath, "supertype_degene_summary_pcut", pthreshold, "_fc", fcthreshold,  "_excludegenes.rds"))

de.results.long = de.results.long0[de.results.long0$padjust < pthreshold & abs(de.results.long0$coef) > fcthreshold,]
saveRDS(de.results.long, file = paste0(opath, "de.results.long_pcut", pthreshold, "_fc", fcthreshold,  ".rds"))











