## Description/Goal:
## Compute UMAP

library(scrattch.bigcat)
library(scrattch.hicat)
library(MatrixGenerics)

#################################################
## Load data
#################################################
load("cluster_0916/merge.result_th300_mincells100.rda") ## clustering results
load("cluster_0916/norm.dat_cluster0916_cl100.rda") ## subsampled cell by gene matrix
load("big.dat_Aging.16ROI.V3.20220701.rda") ## fbm big.dat
load("anno.df.clean_cluster0916_v1020_20221020.rda") ## metadata

## Markers
cl = merge.result$cl
markers = merge.result$markers
select.markers = intersect(markers, row.names(norm.dat))
load("genotype.genes.rda")
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/10X_analysis/sex.genes.rda")
select.markers = setdiff(select.markers, sex.genes)
select.markers = setdiff(select.markers, genotype.genes)


## Rm.eigen
rm.eigen = as.matrix(setNames(log2(anno.df.clean$gene.counts.0), anno.df.clean$sample_id), ncol=1)
rownames(rm.eigen) = anno.df.clean$sample_id
colnames(rm.eigen) = "log2Genes"
save(rm.eigen, file = "rm.eigen.rda") ## save for later use

#################################################
## PCA
#################################################
set.seed(123)

## PCA all
select.cells = anno.df.clean$sample_id
rd.result = rd_PCA_big(big.dat, norm.dat[select.markers, ], select.cells = select.cells, max.dim=150, method="elbow", mc.cores=20)
rd.dat = filter_RD(rd.result$rd.dat, rm.eigen, 0.7)
save(rd.dat, file = "./cluster_0916//rd.result_refdat_2022-10-21.rda") ## save PCA results as R object
write.csv(rd.dat, file = "./cluster_0916/rd.result_refdat_2022-10-21.csv") ## save PCA results as csv which will be used for UMAP computation in Python


#####################################
## UMAP
#####################################
## These lines of script call run_umap.py script for UMAP computation
tmp.in = "./cluster_0916/rd.result_refdat_2022-10-21.csv"
tmp.out = "./cluster_0916/umap.2d_global_k10_d0.4_2022-10-21.csv"
tmp.system = paste0("/functions/run_umap.py -i ", tmp.in, " -o ", tmp.out, " -k 10 -c 2 -d 0.4")
system(tmp.system)