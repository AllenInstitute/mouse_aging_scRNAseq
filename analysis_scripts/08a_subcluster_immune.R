## Description/Goal:
## Extract only immune cells and cluster them separately


library(dplyr)
library(scrattch.hicat)
library(scrattch.bigcat)

#################################################
## Load data
#################################################
load("anno.df.clean_cluster0916_v1020_20221020.rda")
load("big.dat_Aging.16ROI.V3.20220701.rda") ## fbm big.dat

## Create smaller cell-by-gene matrix with only immune cells
select.cells = anno.df.clean$sample_id[anno.df.clean$class_label == "Immune"]
norm.dat = get_logNormal(big.dat, select.cells)



#################################################
## Clustering
#################################################
set.seed(123)

## Cluster
resulti = iter_clust(norm.dat = norm.dat[,select.cells], select.cells = select.cells, "./data_immune/",
                     max.cl.size = 200, split.size = 50, verbose = 1)

## Merge
de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=50, padj.th = 0.01, min.genes = 3) ## Looser settings
merge.result = merge_cl(norm.dat=norm.dat, cl=resulti$cl, rd.dat.t=norm.dat[resulti$markers,select.cells], de.param=de.param, verbose=TRUE) ## 10

## Save for later
save(resulti, file = "./data_immune/resulti.rda")
save(merge.result, file = "./data_immune/merge.result.rda")


