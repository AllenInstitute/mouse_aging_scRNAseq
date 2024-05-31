## Description/Goal:
## Annotate and filter based on latest clustering

library(plyr)
library(dplyr)
library(scrattch.hicat)
options(stringsAsFactors = F)


#################################################
## Load data
#################################################
## Load annotations that already exist from ABC-MWB
load("r_objects/ABC_MWB_taxonomy/AIT21.0_mouse/cl.clean.rda")
remove(cl.df.clean)
load("mapping_aged_fmap_AIT21.0_allcells/ref.cl.df.rda")

## Load clustering results
load("./cluster_0916/merge.result_th300_mincells100.rda")
load("samp.dat_bothages_qcscore_20231208.rda")

cl = merge.result$cl
samp.dat.filtered = samp.dat[samp.dat$sample_id %in% names(cl),]



#################################################
## Combine adult annotations with  cluster results
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

ref.cl.df$cl = as.character(ref.cl.df$cl)
anno.df.u19 = left_join(anno.df.u19, ref.cl.df[,c("cl", "cluster_id_label", "supertype_label", "subclass_label", "class_label")], by = c("cl.u19" = "cl"))
anno.df.u19[is.na(anno.df.u19)] = "missing.or.junk"

head(anno.df.u19)
tail(anno.df.u19)
dim(anno.df.u19)


## Add mapped aged labels
load("mapping_aged_fmap_AIT21.0_allcells//ref.cl.df.rda")
map.result = readRDS("mapping_aged_fmap_AIT21.0_allcells/map.result_bestmap.allcells.20231219.rds")
map.result =  left_join(map.result, ref.cl.df[,c("cl", "cluster_id_label", "supertype_label", "subclass_label", "class_label")], by = c("best.cl_AIT21.0" = "cl"))
names(map.result)[names(map.result) == "best.cl._AIT21.0"] = "cl.u19"
common = intersect(names(map.result), names(anno.df.u19))
map.result = map.result[map.result$sample_id %in% samp.dat.filtered$sample_id[samp.dat.filtered$age_cat == "aged"],]

anno.df.u19 = rbind(anno.df.u19[,common], map.result[,common])


## Add cl to my data
row.names(samp.dat.filtered) = samp.dat.filtered$sample_id
samp.dat.filtered$cl = cl[samp.dat.filtered$sample_id]
samp.dat.filtered = samp.dat.filtered[!is.na(samp.dat.filtered$cl),]

## Add U19 labels
samp.dat.filtered = left_join(samp.dat.filtered, anno.df.u19)

cl.df  = samp.dat.filtered %>% dplyr::group_by(cl) %>% dplyr::summarize(ncell = n(), 
                                                                        med.gc = median(gene.counts.0),
                                                                        med.qc = median(qc.score),
                                                                        med.doub = median(doublet_score),
                                                                        age_prop = sum(age_cat == "aged") / n(),
                                                                        adult_prop = sum(age_cat == "adult") / n(),
                                                                        
                                                                        max.roi = names(which.max(table(roi))),
                                                                        max.roi.prop = (max(table(roi)) / n()),
                                                                        
                                                                        max.class = names(which.max(table(class_label))),
                                                                        max.class_2 = names(sort(table(class_label), decreasing = T))[2],
                                                                        max.class.prop = (max(table(class_label)) / n()),
                                                                        max.class.prop_2 = sort((table(class_label)) / n(), decreasing = T)[2],
                                                                        
                                                                        max.subclass = names(which.max(table(subclass_label))),
                                                                        max.subclass_2 = names(sort(table(subclass_label), decreasing = T))[2],
                                                                        max.subclass.prop = (max(table(subclass_label)) / n()),
                                                                        max.subclass.prop_2 = sort((table(subclass_label)) / n(), decreasing = T)[2],
                                                                        
                                                                        max.supertype = names(which.max(table(supertype_label))),
                                                                        max.supertype_2 = names(sort(table(supertype_label), decreasing = T))[2],
                                                                        max.supertype.prop = (max(table(supertype_label)) / n()),
                                                                        max.supertype.prop_2 = sort((table(supertype_label)) / n(), decreasing = T)[2],
                                                                        
                                                                        max.cluster = names(which.max(table(cluster_id_label))),
                                                                        max.cluster_2 = names(sort(table(cluster_id_label), decreasing = T))[2],
                                                                        max.cluster.prop = (max(table(cluster_id_label)) / n()),
                                                                        max.cluster.prop_2 = sort((table(cluster_id_label)) / n(), decreasing = T)[2],
                                                                        
                                                                        max.lib.prop = (max(table(library_prep)) / n()))



hist(cl.df$max.class.prop)
hist(cl.df$max.subclass.prop)


#################################################
## Vis QC metrics
#################################################
library(ggplot2)
library(gridExtra)


## Gene counts
p2 = cl.df %>%
  ggplot(aes(x = med.gc))+
  geom_histogram(bins = 50)+
  scale_x_log10()+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 3000, color = "red")+
  geom_vline(xintercept = 5000, color = "blue")+
  geom_vline(xintercept = 5500, color = "dark green")+
  ggtitle("Median Gene Counts")+
  facet_wrap(~max.class, ncol = 2, scales = "free_y")+
  theme_bw()

p1 = cl.df %>%
  ggplot(aes(x = med.gc))+
  geom_histogram(bins = 50)+
  scale_x_log10()+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 3000, color = "red")+
  geom_vline(xintercept = 5000, color = "blue")+
  geom_vline(xintercept = 5500, color = "dark green")+
  ggtitle("Median Gene Counts")+
  facet_wrap(~max.class, ncol = 2)+
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
  facet_wrap(~max.class, ncol = 2, scales = "free_y")+
  theme_bw()

p1 = cl.df %>%
  ggplot(aes(x = med.qc))+
  geom_histogram(bins = 50)+
  # scale_x_log10()+
  geom_vline(xintercept = 150)+
  geom_vline(xintercept = 300, color = "red")+
  ggtitle("Median QC Score")+
  facet_wrap(~max.class, ncol = 2)+
  theme_bw()

grid.arrange(p1, p2, nrow = 1)




#################################################
## Label LQ clusters
#################################################
## Cluster-level QC and GC cutoffs for different classes
tmp = c("Immune", "Vascular", "Astro-Epen")
tmp.nn = cl.df$cl[cl.df$max.class %in% tmp & (cl.df$med.gc < 2000 | cl.df$med.qc < 150)]
tmp = c("OPC-Oligo")
tmp.oligo = cl.df$cl[cl.df$max.class %in% tmp & (cl.df$med.gc < 3000 | cl.df$med.qc < 150)]
tmp = c("OB-CR Glut", "DG-IMN Glut", "OB-IMN GABA")
tmp.imn = cl.df$cl[cl.df$max.class %in% tmp & (cl.df$med.gc < 3000 | cl.df$med.qc < 150)]
tmp = setdiff(unique(ref.cl.df$class_label), c("Immune", "Vascular", "OPC-Oligo", "Astro-Epen", "OB-CR Glut", "DG-IMN Glut", "OB-IMN GABA"))
tmp.neuron = cl.df$cl[cl.df$max.class %in% tmp & (cl.df$med.gc < 5500 | cl.df$med.qc < 300)]

cl.lq.qc = unique(c(tmp.nn, tmp.imn, tmp.neuron, tmp.oligo))

## Clusters that are mostly mapping to "junk"
cl.lq.junk = cl.df$cl[cl.df$max.subclass == "missing.or.junk" | cl.df$max.supertype == "missing.or.junk"]

## Clusters that are biased towards a single donor
cl.lq.donor = cl.df$cl[cl.df$max.donor.prop > 0.8]

## Other manually annotated LQ clusters:
cl.lq.manual = c("2108" #mapped to half neurons and half microglia
)

## Exclude off target supertypes
load("r_objects/supertype.exclude.rda")
cl.offtarget = cl.df$cl[cl.df$max.supertype %in% supertype.exclude]
tmp = cl.df[cl.df$cl %in% cl.offtarget,]

## Ambiguous mapping
hist(cl.df$max.class.prop)
tmp.cl.df = cl.df[cl.df$max.class.prop < 0.7,]
cl.ambmappingjunk = tmp.cl.df$cl[tmp.cl.df$max.class_2 == "missing.or.junk" & tmp.cl.df$max.class.prop_2 > 0.3]
tmp = cl.df[cl.df$cl %in% unique(c(cl.ambmappingjunk)),]

cl.lq = unique(c(cl.lq.qc, cl.lq.junk, cl.lq.donor, cl.lq.manual, cl.offtarget, cl.ambmappingjunk))


#################################################
## Final filter
#################################################
cl.df$keep = T
cl.df$keep[cl.df$cl %in% cl.lq] = F

sum(cl.df$keep) ## number of remaining clusters
sum(cl.df$ncell) ## number of starting cells
sum(cl.df$ncell[cl.df$keep == T]) ## number of cells remaining after new filters

cl.df.clean = cl.df[cl.df$keep == T,]

save(cl.df, file = "./cluster_0916/cl.df_v231221.rda")

## Assign class
cl.df.clean$class_label = cl.df.clean$max.class
table(cl.df.clean$class_label)

## Assign subclass
cl.df.clean$subclass_label = cl.df.clean$max.subclass
table(cl.df.clean$subclass_label)

## Assign supertype
cl.df.clean$supertype_label = cl.df.clean$max.supertype
table(cl.df.clean$supertype_label)


#################################################
## Finding hierarchy discrepencies - correct each manually
#################################################
load("mapping_aged_fmap_AIT21.0_allcells//ref.cl.df.rda")
tmp = unique(cl.df.clean[,c("class_label", "subclass_label", "supertype_label")])
tmp.ref = unique(ref.cl.df[,c("class_label", "subclass_label", "supertype_label")])
setdiff(tmp, tmp.ref)

##--> cl292
cl.df.clean$class_label[cl.df.clean$cl == 292] = "MB Glut" #changing to MB Glut because its similar in prop to former (P Glut) and max ROI is also MB
##--> cl502
cl.df.clean$supertype_label[cl.df.clean$cl == 502] = "RHP-COA Ndnf Gaba_6" #changing to top supertype of the matching subclass
##--> cl629
cl.df.clean$supertype_label[cl.df.clean$cl == 629] = "STN-PSTN Pitx2 Glut_2" #changing to top supertype of the matching subclass
##--> cl3678
cl.df.clean$supertype_label[cl.df.clean$cl == 3678] = "L4/5 IT CTX Glut_1" #changing to top supertype of the matching subclass
##--> cl3415
cl.df.clean$supertype_label[cl.df.clean$cl == 3415] = "OB-STR-CTX Inh IMN_1" #changing to top supertype of the matching subclass; all aged cells in this cluster are also mapping to OB rather than DIR

tmp = unique(cl.df.clean[,c("class_label", "subclass_label", "supertype_label")])
tmp.ref = unique(ref.cl.df[,c("class_label", "subclass_label", "supertype_label")])
setdiff(tmp, tmp.ref) #should be empty now

## Correcting order
class_label = intersect(unique(ref.cl.df$class_label), unique(cl.df.clean$class_label))
class_id = 1:length(class_label)
class.df = data.frame(class_label = class_label, class_id = class_id)

subclass_label = intersect(unique(ref.cl.df$subclass_label), unique(cl.df.clean$subclass_label))
subclass_id = 1:length(subclass_label)
subclass.df = data.frame(subclass_label = subclass_label, subclass_id = subclass_id)

supertype_label = intersect(unique(ref.cl.df$supertype_label), unique(cl.df.clean$supertype_label))
supertype_id = 1:length(supertype_label)
supertype.df = data.frame(supertype_label = supertype_label, supertype_id = supertype_id)

cl.df.clean = left_join(cl.df.clean, class.df)
cl.df.clean = left_join(cl.df.clean, subclass.df)
cl.df.clean = left_join(cl.df.clean, supertype.df)


## Add colors
class.colors = read.csv("mapping_aged_fmap_AIT21.0/ait21_class_colors.csv")
subclass.colors = read.csv("mapping_aged_fmap_AIT21.0/ait21_subclass_colors.csv")
supertype.colors = read.csv("mapping_aged_fmap_AIT21.0/ait21_supertype_colors.csv")

cl.df.clean = left_join(cl.df.clean, class.colors)
cl.df.clean = left_join(cl.df.clean, subclass.colors)
cl.df.clean = left_join(cl.df.clean, supertype.colors)


## Update cl numbers
cl.df.clean = cl.df.clean[order(cl.df.clean$supertype_id, cl.df.clean$subclass_id, cl.df.clean$class_id),]
cl.df.clean$cl_0916 = cl.df.clean$cl
cl.df.clean$cl = 1:nrow(cl.df.clean)
cl.df.clean$cluster_label = paste0(cl.df.clean$cl, "_", cl.df.clean$supertype_label)

save(cl.df.clean, file = "cluster_0916/cl.df.clean_v231221.rda")


#################################################
## Make new annotated metadata object (anno.df.clean)
#################################################
## Load clustering results and new cl.df.clean
load("./cluster_0916/merge.result_th300_mincells100.rda")
load("./cluster_0916/cl.df.clean_v231221.rda")

## Load original cell metadata file and key for broad rois
load("samp.dat_bothages_qcscore_20231208.rda")
load("r_objects/broad.roi.key.rda")

## Add newest cluster label to metadata file
cl = merge.result$cl
select.cells = names(cl)
row.names(samp.dat) = samp.dat$sample_id
samp.dat.filtered = samp.dat[select.cells,]
samp.dat.filtered$cl_0916 = cl[row.names(samp.dat.filtered)]

anno.df.clean = samp.dat.filtered[samp.dat.filtered$cl_0916 %in% cl.df.clean$cl_0916,]
anno.df.clean = left_join(anno.df.clean, cl.df.clean)
anno.df.clean = left_join(anno.df.clean, broad.key)

save(anno.df.clean, file = "anno.df.clean_cluster0916_v231221_20231221.rda") 

