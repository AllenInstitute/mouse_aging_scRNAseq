## Description/Goal:
## Running MAST at the subclass level to find age-DE genes

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

## Scrattch.bigcat compatible file-backed-matrix
load("big.dat_Aging.16ROI.V3.20231208.rda")



################################################################################
## Deciding which subclasses to include for analysis
# ################################################################################
# -- Criteria:
# -- ncell > 100 for each age_cat

tmp.age.count = c(table(anno.df.clean$subclass_label[anno.df.clean$age_cat == "aged"]))
tmp.adult.count = c(table(anno.df.clean$subclass_label[anno.df.clean$age_cat == "adult"]))
tmp.age.count2 = tmp.age.count[tmp.age.count > 100]
tmp.adult.count2 = tmp.adult.count[tmp.adult.count > 100]

tmp.keep = intersect(names(tmp.age.count2), names(tmp.adult.count2))

check.remove = setdiff(unique(anno.df.clean$subclass_label), tmp.keep)

subclass.keep = tmp.keep

save(subclass.keep, file = "subclass.keep_cluster0916_v20231221.rda")



################################################################################
## MAST per subclass (using parallel function)
################################################################################
source("functions/degenes_age_mast_simplified_qcscore.R") ## custom function for running MAST in parallel
setwd("MAST_simplified/")


de.gene.results = get_age_degenes_mast_big_p(anno.df.clean, big.dat, groups = subclass.keep, split.by = "subclass_label", save.as.tmp = F,
                                    maxcells = 2000, freq_expressed = 0.1, mc.cores = 15)

saveRDS(de.gene.results, file = "de.gene.result_cluster0916_v231221_subclass.rds")



################################################################################
## Assemble MAST results
################################################################################
## Load results
load("anno.df.clean_cluster0916_v231221_20240206.rda")
results = readRDS("MAST_simplified/de.gene.result_cluster0916_v231221_subclass.rds")

## Exclude genes
load("r_objects/sex.genes.rda") #sex-biased genes
load("r_objects/genotype.genes.rda") #genotype-biased genes
load("r_objects/iegenes.rda") #immediate early genes
genes.exclude = c(iegenes.rapidprg, iegenes.delayedprg, sex.genes, genotype.genes)

## Cutoffs for age-DE genes
fcthreshold = 1
pthreshold = 0.01


## Create summarized objects with results
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

saveRDS(subclass.results, file = paste0(opath, "subclasswide_degene_summary_pcut", pthreshold, "_fc", fcthreshold,  "_excludedgenes.rds"))

de.results.long = de.results.long0[de.results.long0$padjust < pthreshold & abs(de.results.long0$coef) > fcthreshold,]
saveRDS(de.results.long, file = paste0(opath, "de.results.long_pcut", pthreshold, "_fc", fcthreshold,  ".rds"))


################################################################################
## Visualize results
################################################################################
opath = "MAST_simplified"

tmp.level1 = unique(anno.df.clean[,c("class_label", "class_id")])
tmp.level1 = tmp.level1[order(tmp.level1$class_id),]
level1.order = tmp.level1$class_label
tmp = unique(anno.df.clean[,c("subclass_label", "class_color", "class_label")])
level1.color = tmp$class_color
names(level1.color) = tmp$subclass_label

subclass.results1 = left_join(subclass.results, tmp, by = c("subclass" = "subclass_label"))
row.names(subclass.results1) = subclass.results1$subclass
# subclass.results1$subclass_color = level1.color[subclass.results1$subclass]
subclass.results1 = subclass.results1[order(subclass.results1$subclass),]
subclass.results1$class_label = factor(subclass.results1$class_label, levels = level1.order)
subclass.results1 = subclass.results1[order(subclass.results1$class_label),]
subclass.results1$subclass = factor(subclass.results1$subclass, levels = subclass.results1$subclass)

## Barplot of nDEgene count
p1 = subclass.results1 %>%
  ggplot(aes(y = de.gene.count, x = subclass, fill = subclass))+
  geom_bar(stat = "identity", width = 0.8)+
  scale_fill_manual(values = level1.color)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))


de.results.long$subclass = factor(de.results.long$subclass, levels = subclass.results1$subclass)
de.results.long$direction = "positive"
de.results.long$direction[de.results.long$coef < 0] = "negative"
de.results.long$abs.coef = abs(de.results.long$coef)


## Dot plot of coeficients
p2 = de.results.long %>%
  ggplot(aes(x = subclass, y = abs.coef, color = subclass, shape = direction))+
  geom_jitter(size = 0.8, width = 0.05, alpha = 0.8)+
  scale_color_manual(values = level1.color)+
  scale_shape_manual(values = c(4, 1))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

p4 = de.results.long %>%
  ggplot(aes(x = subclass, y = -log10(padjust), color = subclass, size = 5))+
  geom_jitter(size = 0.5, width = 0.1)+
  scale_color_manual(values = level1.color)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))


p3 = subclass.results1 %>%
  ggplot(aes(y = de.gene.count, x = subclass, fill = subclass))+
  geom_bar(stat = "identity", width = 0.8)+
  scale_fill_manual(values = level1.color)+
  geom_jitter(data = de.results.long, aes(x = subclass, y = -30*abs.coef, color = subclass, shape = direction), size = 0.8, width = 0.05, alpha = 0.8)+
  scale_shape_manual(values = c(4, 1))+
  scale_color_manual(values = level1.color)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))


ggsave(p1, filename = paste0(opath, "barplot_ndegene_level1color_pcut", pthreshold, "_fc", fcthreshold, ".pdf"), useDingbats = F, width = 14, height = 5)
ggsave(p2, filename = paste0(opath, "dotplot_signifcoef_level1color_pcut", pthreshold, "_fc", fcthreshold, ".pdf"), useDingbats = F, width = 14, height = 5)
ggsave(p4, filename = paste0(opath, "dotplot_signifpval_level1color_pcut", pthreshold, "_fc", fcthreshold, ".pdf"), useDingbats = F, width = 14, height = 5)
ggsave(p3, filename = paste0(opath, "combdotbar_level1color_pcut", pthreshold, "_fc", fcthreshold, ".pdf"), useDingbats = F, width = 14, height = 5)
