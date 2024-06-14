## Description/Goal:
## Run pseudobulk analysis on all subclasses


###################
## Load data
###################
library(scrattch.hicat)
library(scrattch.bigcat)
library(dplyr)

load("anno.df.clean_cluster0916_v230915_20230915.rda")
load("cluster_0916/norm.dat_cluster0916_v230915_subclass_ss1500.rda") ## subsampled cell by gene matrix
load("big.dat_Aging.16ROI.V3.20220701.rda")
load("subclass.keep_cluster0916_v230915.rda")

###################
## Compute pseudobulk means (means per library per subclass)
###################
anno.df.select = anno.df.clean[anno.df.clean$subclass_label %in% subclass.keep,]
cl = paste(anno.df.select$subclass_label, anno.df.select$age_cat, anno.df.select$library_prep, sep = "__")
length(unique(cl))
names(cl) = anno.df.select$sample_id

tmp.summary = anno.df.select %>% dplyr::group_by(subclass_label, age_cat, library_prep) %>% dplyr::summarize(n = n())
tmp.summary.keep = tmp.summary[tmp.summary$n >= 10,]
tmp.summary.keep$cl = paste(tmp.summary.keep$subclass_label, tmp.summary.keep$age_cat, tmp.summary.keep$library_prep, sep = "__")

cl.clean = cl[cl %in% tmp.summary.keep$cl]

cl.means = get_cl_stats_big(big.dat, cl.clean, max.cl.size = 200, stats = "means", mc.cores = 5)

save(cl.means, file = "pb_edgeR/cl.means")



###################
## Pseudobulk analysis with EdgeR
###################
library(edgeR)
library(statmod)

load("pb_edgeR/cl.means")
load("subclass.keep_cluster0916_v230915.rda")

cl.means = cl.means$means
tmp.split = strsplit(colnames(cl.means), "__")
subclass_label = unlist(lapply(tmp.split, function(x) x[[1]]))
age_cat = unlist(lapply(tmp.split, function(x) x[[2]]))
library_prep = unlist(lapply(tmp.split, function(x) x[[3]]))
meta.df = data.frame(sample = colnames(cl.means), subclass_label = subclass_label, age_cat = age_cat, library_prep = library_prep)

qlf.result = c()
lrt.result = c()

for(i in subclass.keep){
  
  print(i)
  
  sample.keep = grep(i, colnames(cl.means))
  meta.select = meta.df[sample.keep,]
  
  age = meta.select$age_cat
  # y.select = DGEList(counts = cl.means[,sample.keep], group = age)
  
  tmp.counts = 2^cl.means[,sample.keep] - 1
  y.select.counts = DGEList(counts = tmp.counts, group = age)
  
  design = model.matrix(~age)
  row.names(design) = colnames(y.select.counts)
  
  y.select.counts = estimateDisp(y.select.counts, design, robust = T)
  # plotBCV(y.select.counts)
  
  print("Running quasi-likelihood...")
  fit = glmQLFit(y.select.counts, design)
  qlf = glmQLFTest(fit, coef = 2)
  # topTags(qlf, n = 20)
  
  print("Running LRT...")
  fit = glmFit(y.select.counts, design)
  lrt = glmLRT(fit, coef = 2)
  # topTags(lrt, n = 20)
  
  print("Saving results...")
  tmp.qlf = list(result = qlf, dgelist = y.select.counts)
  tmp.lrt = list(result = lrt, dgelist = y.select.counts)
  qlf.result[i] = list(tmp.qlf)
  lrt.result[i] = list(tmp.lrt)
  
}


saveRDS(qlf.result, file = "pb_edgeR/qlf.result_subclass_cluster0916_v230915.rds")
saveRDS(lrt.result, file = "pb_edgeR/lrt.result_subclass_cluster0916_v230915.rds")



