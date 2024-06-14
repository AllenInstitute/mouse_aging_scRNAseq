## Description/Goal:
## Gene set enrichment analysis with gprofiler2 at the subclass level

library(gprofiler2)
library(dplyr)

#################################################
## Run gprofiler
#################################################
## Load in age-DE genes from MAST
de.results = readRDS("MAST_simplified/de.gene.result_cluster0916_v231221_subclass.rds")
de.results = lapply(de.results, function(x) x[["result"]])

## Remove certain genes
load("sex.genes.rda")
load("genotype.genes.rda")
load("iegenes.rda")
genes.exclude = c(iegenes.rapidprg, iegenes.delayedprg, sex.genes, genotype.genes)

de.results = lapply(de.results, function(x) x[!x$primerid %in% genes.exclude,])

## Signif cutoffs
pval.cut = 0.01
logfc.cut = 1

## Positive only
de.genes = lapply(de.results, function(x) x[abs(x$coef) > logfc.cut & x$padjust < pval.cut & x$coef > 0,])
tmp = unlist(lapply(de.genes, nrow))
de.genes = de.genes[names(tmp[tmp >= 3])]
gp.result.pos = lapply(de.genes, function(x) gost(query = x$primerid, organism = "mmusculus", ordered_query = T, evcodes = T, user_threshold = 0.01))


## Negative only
de.genes = lapply(de.results, function(x) x[abs(x$coef) > logfc.cut & x$padjust < pval.cut & x$coef < 0,])
tmp = unlist(lapply(de.genes, nrow))
de.genes = de.genes[names(tmp[tmp >= 3])]
gp.result.neg = lapply(de.genes, function(x) gost(query = x$primerid, organism = "mmusculus", ordered_query = T, evcodes = T, user_threshold = 0.01))


save(gp.result.pos, gp.result.neg, file = "gprofiler_results_subclass_all_AIT21_logFC1_pval0.01cut_posneg_separate_excludedgenes_20240109.rda")



#################################################
## Consolidate results
#################################################
load("gprofiler_results_subclass_all_AIT21_logFC1_pval0.01cut_posneg_separate_excludedgenes_20240109.rda")

for(i in names(gp.result.pos)){
  print(i)
  tmp = gp.result.pos[[i]]$result
  if(is.null(tmp)){
    gp.result.pos[[i]] = NULL
    next
  }
  tmp$subclass = i
  tmp$direction = "pos"
  gp.result.pos[[i]] = tmp
}

for(i in names(gp.result.neg)){
  print(i)
  tmp = gp.result.neg[[i]]$result
  if(is.null(tmp)){
    gp.result.neg[[i]] = NULL
    next
  }
  tmp$subclass = i
  tmp$direction = "neg"
  gp.result.neg[[i]] = tmp
}

tmp.pos = do.call(rbind, gp.result.pos)
tmp.neg = do.call(rbind, gp.result.neg)

go.results.all = rbind(data.frame(tmp.pos), data.frame(tmp.neg))
go.results.all = apply(go.results.all,2,as.character)

write.csv(go.results.all, file = "gprofiler_results_subclass_all_AIT21_logFC1_pval0.01cut_posneg_separate_excludedgenes_20240109.csv")