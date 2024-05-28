## Description/Goal:
## Perform final round of clustering


library(dplyr)
library(scrattch.hicat)
library(scrattch.bigcat)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)


#################################################
## Load data
#################################################
load("big.dat_Aging.16ROI.V3.20220701.rda")
load("anno.df.clean_20220915.rda") ## most recent metadata object with filtered cells

select.cells = anno.df.clean$sample_id


#################################################
## Second clustering round
#################################################
de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=50)

set.seed(123)
resulti = iter_clust_big(big.dat = big.dat, select.cells = select.cells, de.param =de.param, prefix="./cluster_0916/cluster_intermed_files",
                         max.cl.size=200, split.size = 50, verbose=1, sampleSize=10000)

save(resulti, file=paste0("./cluster_0916/resulti_cluster.rda"))


## Merge clusters
load("./cluster_0916/resulti_cluster.rda")
cl = resulti$cl
markers = resulti$markers

sampled.cells = sample_cells(cl, 100)
# sampled.cells = sample(sampled.cells, 200000, replace = F)
norm.dat = get_logNormal(big.dat, sampled.cells)
select.markers = intersect(markers, big.dat$row_id)

de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=100)
merge.result = merge_cl(norm.dat=norm.dat, cl=cl, rd.dat.t=norm.dat[select.markers,], de.param=de.param, verbose=TRUE)



#################################################
## Save objects for future use
#################################################
save(merge.result, file= "./cluster_0916/merge.result_th300_mincells100.rda") ## clustering results
save(norm.dat, file = "./cluster_0916/norm.dat_cluster0916_cl100.rda") ## subsampled cell-by-gene matrix
