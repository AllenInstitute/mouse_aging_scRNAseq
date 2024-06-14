## Description/Goal:
## Run augur on subclasses

library(Augur)

#################################################
## Modify select_variance function in Augur
#################################################
library(godmode) ## this package used for running modified version of select_variance()
source("cust_functions/select_variance_MOD.R")
godmode:::assignAnywhere("select_variance", select_variance2)

#################################################
## Load data and combine
#################################################
load("anno.df.clean_cluster0916_v231221_20231221.rda")
load("subclass.keep_cluster0916_v20231221.rda") ## subclasses included in MAST analysis
load("big.dat_Aging.16ROI.V3.20231208.rda")
load("cluster_0916/norm.dat_cluster0916_v20231221_subclass_ss1200.rda")

#################################################
## Run Augur
#################################################
anno.df.select = anno.df.clean[anno.df.clean$sample_id %in% colnames(norm.dat),]
anno.df.select = anno.df.select[anno.df.select$subclass_label %in% subclass.keep,]
anno.df.select$subclass_label = as.character(anno.df.select$subclass_label)
min(table(anno.df.select$subclass_label, anno.df.select$age_cat))

meta.df = anno.df.select[,c("subclass_label", "age_cat")]
row.names(meta.df) = anno.df.select$sample_id

norm.dat.select = norm.dat[,row.names(meta.df)]

augur = calculate_auc(norm.dat.select, meta.df, cell_type_col = "subclass_label", label_col = "age_cat", 
                      feature_perc = 0.8, var_quantile = 0.9, subsample_size = 200,
                      n_threads = 25)

saveRDS(augur, file = "augur/augur.result_ait21_dispersioncorrected_var0.9_fperc0.8_subsampsize200_20240102.rds")

#################################################
## Check results and save smaller objects for later
#################################################
augur = readRDS("augur/augur.result_ait21_dispersioncorrected_var0.9_fperc0.8_subsampsize200_20240102.rds")
rank.all = augur$AUC
aucresults.all = augur$results
save(rank.all, aucresults.all, file = "augur/augur.objects_allcells.rda")

load("augur/augur.objects_allcells.rda")

rank.all = rank.all[order(rank.all$cell_type),]
names(rank.all) = c("cell_type", "auc.all")
rank.summary = rank.summary[order(rank.summary$auc.all, decreasing = T),]
saveRDS(rank.summary,  file = "augur/rank.summary_all.m.f_ait21_dispersioncorrected_var0.9_fperc0.8_subsampsize200_20240102.rds")
