## Description/Goal:
## Use preliminary clustering results to perform additional round of filtering of clusters 
## before second and final round of clustering

## Note:
## This script uses an older version of ABC-MWB annotations (AIT13.0) - we recommend using
## the latest version that has been published in Yao, 2023 Nature, but are including the older
## version here for transparancy


library(dplyr)
library(scrattch.hicat)
options(stringsAsFactors = F)


#################################################
## Load data files
#################################################
## Load objects from prior scripts
load("./prelim_cluster_objects/merge.result_th300_cluster_0220705.rda")
load("samp.dat.filtered_prelim.20220705.rda")

## Load annotations that already exist for adult cells (these were an older version of ABC-MWB annotations)
load("r_objects/ABC_MWB_taxonomy/AIT13.0_mouse/cl.clean.rda") ## loads objects cl.clean and cl.df.clean


#################################################
## Combine adult annotations with preliminary cluster results
#################################################
## Subset ABC-MWB results and label missing cells (assumed to be junk)
roi.keep = unique(samp.dat.filtered$roi)

cells.anno = intersect(samp.dat.filtered$sample_id, names(cl.clean))
cells.adult = samp.dat.filtered$sample_id[samp.dat.filtered$age_cat == "adult"]
cells.junk = setdiff(cells.adult, cells.anno)

anno.df.u19 = data.frame(sample_id = c(cells.anno, cells.junk))
anno.df.u19$cl.u19 = cl.clean[anno.df.u19$sample_id]
anno.df.u19$cl.u19[is.na(anno.df.u19$cl.u19)] = "missing.or.junk"

head(anno.df.u19)
tail(anno.df.u19)

anno.df.u19 = left_join(anno.df.u19, cl.df.clean[,c("cl", "cluster_label", "Level1_label", "Level2_label", "subclass_label", "class_label")], by = c("cl.u19" = "cl"))
anno.df.u19$cluster_label = as.character(anno.df.u19$cluster_label)
anno.df.u19[is.na(anno.df.u19)] = "missing.or.junk"

head(anno.df.u19)
tail(anno.df.u19)


## Add cl to my data
row.names(samp.dat.filtered) = samp.dat.filtered$sample_id
cl = merge.result$cl
samp.dat.filtered$cl = cl[samp.dat.filtered$sample_id]

## Add U19 labels
samp.dat.filtered = left_join(samp.dat.filtered, anno.df.u19)


cl.df  = samp.dat.filtered %>% group_by(cl) %>% summarize(ncell = n(), ## this object summarizes many QC stats for each de novo cluster and will be used to visualize and determine filter criteria
                                                          med.gc = median(gene.counts.0),
                                                          med.qc = median(qc.score),
                                                          med.doub = median(doublet_score),
                                                          age_prop = sum(age_cat == "aged") / n(),
                                                          adult_prop = sum(age_cat == "adult") / n(),
                                                          U19_id_prop = 1 - (sum(is.na(cluster_label)) / n()),
                                                          max.cluster.U19 = names(which.max(table(cluster_label))),
                                                          max.cluster.U19.prop = (max(table(cluster_label)) / n()),
                                                          max.subclass.U19 = names(which.max(table(subclass_label))),
                                                          max.subclass.U19.prop = (max(table(subclass_label)) / n()),
                                                          max.class.U19 = names(which.max(table(class_label))),
                                                          max.class.U19.prop = (max(table(class_label)) / n()),
                                                          max.level1.U19 = names(which.max(table(Level1_label))),
                                                          max.level1.U19.prop = (max(table(Level1_label)) / n()),
                                                          max.level2.U19 = names(which.max(table(Level2_label))),
                                                          max.level2.U19.prop = (max(table(Level2_label)) / n()),
                                                          max.donor.prop = (max(table(library_prep)) / n()))


#################################################
## Vis QC metrics
#################################################
library(gridExtra)
library(ggplot2)

## Gene counts
p2 = cl.df %>%
  ggplot(aes(x = med.gc))+
  geom_histogram(bins = 50)+
  scale_x_log10()+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 3000, color = "red")+
  geom_vline(xintercept = 5500, color = "blue")+
  ggtitle("Median Gene Counts")+
  facet_wrap(~max.class.U19, ncol = 2, scales = "free_y")+
  theme_bw()

p1 = cl.df %>%
  ggplot(aes(x = med.gc))+
  geom_histogram(bins = 50)+
  scale_x_log10()+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 3000, color = "red")+
  geom_vline(xintercept = 5500, color = "blue")+
  ggtitle("Median Gene Counts")+
  facet_wrap(~max.class.U19, ncol = 2)+
  theme_bw()

grid.arrange(p1, p2, nrow = 1)

## QC Score
p2 = cl.df %>%
  ggplot(aes(x = med.qc))+
  geom_histogram(bins = 50)+
  # scale_x_log10()+
  geom_vline(xintercept = 150)+
  geom_vline(xintercept = 300, color = "red")+
  ggtitle("Median QC Score")+
  facet_wrap(~max.class.U19, ncol = 2, scales = "free_y")+
  theme_bw()

p1 = cl.df %>%
  ggplot(aes(x = med.qc))+
  geom_histogram(bins = 50)+
  # scale_x_log10()+
  geom_vline(xintercept = 150)+
  geom_vline(xintercept = 300, color = "red")+
  ggtitle("Median QC Score")+
  facet_wrap(~max.class.U19, ncol = 2)+
  theme_bw()

grid.arrange(p1, p2, nrow = 1)


#################################################
## Label LQ clusters
#################################################
## Cluster-level QC and GC cutoffs for neurons versus non-neurons
tmp1 = cl.df$cl[cl.df$max.class.U19 %in% c("NN", "IMN") & (cl.df$med.gc < 2000 | cl.df$med.qc < 100)] 
tmp.neurons = c("Gaba", "Glut", "Hybrid", "Sero", "Dopa")
tmp2 = cl.df$cl[cl.df$max.class.U19 %in% tmp.neurons & (cl.df$med.gc < 3000 | cl.df$med.qc < 250)]
cl.lq.qc = unique(c(tmp1, tmp2))

## Clusters that are mostly mapping to "junk"
cl.lq.junk = cl.df$cl[cl.df$max.cluster.U19 == "missing.or.junk" & cl.df$U19_id_prop > 0.05]

## Clusters that are biased towards a single donor
cl.lq.donor = cl.df$cl[cl.df$max.donor.prop > 0.8]

## Other manually annotated LQ clusters:
cl.lq.manual = c(
                 7051, ## Oligo cluster with endo markers
                 7191, 7193, ## Astro + Microglia doublets
                 7214, 7215, ## Astro + Endo doublets
                 595 ## Astro + Gluta doublet (Slc17a7)
                 )


cl.lq = unique(c(cl.lq.qc, cl.lq.junk, cl.lq.donor, cl.lq.manual))

cl.df$keep = T
cl.df$keep[cl.df$cl %in% cl.lq] = F

sum(cl.df$keep) ## number of remaining clusters
sum(cl.df$ncell) ## number of starting cells
sum(cl.df$ncell[cl.df$keep == T]) ## number of cells remaining after new filters

save(cl.df, file = "./prelim_cluster_objects/cl.df_U19.max.2022-09-15.rda") ## save cl.df object for later


#################################################
## Filter LQ clusters
#################################################
## Filter
cl = merge.result$cl
samp.dat.filtered$cl = cl[samp.dat.filtered$sample_id]
samp.dat.filtered = left_join(samp.dat.filtered, cl.df)

anno.df.clean = samp.dat.filtered[samp.dat.filtered$keep == T,]

save(anno.df.clean, file = "anno.df.clean_20220915.rda")




## Create anno.df (remove weird geno/donors)
weird.donor = "483935"
weird.geno = c("Gad2-IRES-Cre/wt;Ai14(RCL-tdT)/wt", "CX3CR1-GFP/CX3CR1-GFP")
anno.df.clean = samp.dat.filtered[samp.dat.filtered$external_donor_name != weird.donor,]
anno.df.clean = anno.df.clean[!anno.df.clean$full_genotype %in% weird.geno,]
row.names(anno.df.clean) = anno.df.clean$sample_id
save(anno.df.clean, file = "anno.df.clean_20220915.rda")
