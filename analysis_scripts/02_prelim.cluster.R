## Description/Goal:
## Perform first round of clustering


library(dplyr)
library(scrattch.hicat)
library(scrattch.bigcat)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)


#################################################
## Load data files
#################################################
load("big.dat_Aging.16ROI.V3.20220701.rda") ## cell by gene file-backed matrix (see github.com/alleninstitute/scrattch.bigcat for examples for assembly)
load("samp.dat.filtered_prelim.20220705.rda") ## filtered metadata file from previous step

select.cells = samp.dat.filtered$sample_id


#################################################
## First clustering round
#################################################
## Run clustering
de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=50) ## Looser settings

set.seed(123)
resulti = iter_clust_big(big.dat = big.dat, select.cells = select.cells, de.param =de.param, prefix="./prelim_cluster/cluster_0220705",
                         max.cl.size=200, split.size = 50, verbose=1, sampleSize=10000)

save(resulti, file="./prelim_cluster/resulti_cluster_0220705.rda")

cl = resulti$cl ## de novo clusters
markers = resulti$markers ## marker genes for all clusters

## Subsample results and merge
sampled.cells = sample_cells(cl, 100)
sampled.cells = sample(sampled.cells, 200000, replace = F)
norm.dat = get_logNormal(big.dat, sampled.cells)
select.markers = intersect(markers, big.dat$row_id)

de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=50)
merge.result = merge_cl(norm.dat=norm.dat, cl=cl, rd.dat.t=norm.dat[select.markers,], de.param=de.param, verbose=TRUE)



#################################################
## Save objects for future use
#################################################
save(merge.result, file="merge.result_th300_cluster_0220705.rda")
save(norm.dat, file = "norm.dat_cl100_ss200Kcells.rda")
