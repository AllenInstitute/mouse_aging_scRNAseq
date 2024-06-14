## Description/Goal:
## Compute and plot gene enrichment matrix

library(pheatmap)
library(RColorBrewer)

#################################################
## Load data and combine
#################################################
go.results.all = read.csv("gprofiler_results_subclass_all_AIT21_logFC1_pval0.01cut_posneg_separate_excludedgenes_20240109.csv")
go.results.nn = read.csv("gprofiler_results_nnsupertype_AIT21_logFC1_pval0.01cut_posneg_separate_excludedgenes_20240109.csv")
go.results.all = go.results.all[!grepl(" NN", go.results.all$subclass),]
go.results.all0 = rbind(go.results.all, go.results.nn)


#################################################
## Filter
#################################################
## Which data sources to use
source.keep = c("GO:MF", "GO:BP", "GO:CC")
go.results.all = go.results.all0[go.results.all0$source %in% source.keep,]

## Which terms/subclasses/supertypes to keep
termfreq = 2 #filters terms with frequency less than this
subclassfreq = 6 #filters celltypes with frequency less than this
goterm.keep = names(table(go.results.all$term_name)[table(go.results.all$term_name) >= termfreq])
go.results.all = go.results.all[go.results.all$term_name %in% goterm.keep,]
subclass.keep = names(table(go.results.all$subclass)[table(go.results.all$subclass) >= subclassfreq])
go.results.all = go.results.all[go.results.all$subclass %in% subclass.keep,]
setdiff(unique(go.results.all0$subclass), subclass.keep) # check subclasses that got filtered out

## Filter out broad terms
nterm_cutoff = 2000 #filters terms with greater than this number of gene members (to filter out terms that are too broad)
terms.tmp = unique(go.results.all[,c("term_name", "term_size")])
unique(go.results.all$term_name[go.results.all$term_size > nterm_cutoff]) #which terms would be thrown out
go.results.all = go.results.all[go.results.all$term_size <= nterm_cutoff,]

## Final counts of terms/subclasses
length(unique(go.results.all$term_name))
length(unique(go.results.all$subclass))

go.mat = matrix(0, nrow = length(unique(go.results.all$subclass)), ncol = length(unique(go.results.all$term_name)))
row.names(go.mat) = unique(go.results.all$subclass)
colnames(go.mat) = unique(go.results.all$term_name)

for(n in 1:nrow(go.results.all)){
  x = go.results.all$subclass[n]
  y = go.results.all$term_name[n]
  tmp = -log10(go.results.all$p_value[n])
  if(go.results.all$direction[n] == "pos"){
    go.mat[x,y] = tmp
  }
  if(go.results.all$direction[n] == "neg"){
    go.mat[x,y] = -tmp
  }
}


#################################################
## Plot
#################################################
## Make objects for colors
load("col.broi.rda")
tmp = "#808080"
names(tmp) = "NN"
col.broi = c(col.broi, tmp)
subclass.roi.summary = anno.df.clean %>% group_by(subclass_label) %>% summarize(max.broad.roi = names(which.max(table(broad_roi))))
subclass.roi.summary = subclass.roi.summary[!grepl("NN", subclass.roi.summary$subclass_label),]
tmp = data.frame(subclass_label = unique(go.results.nn$subclass), max.broad.roi = "NN")
subclass.roi.summary = rbind(subclass.roi.summary, tmp)
subclass.roi.summary = subclass.roi.summary[subclass.roi.summary$subclass_label %in% subclass.keep,]

anno.df = data.frame(max.broad.roi = subclass.roi.summary$max.broad.roi)
row.names(anno.df) = subclass.roi.summary$subclass_label

anno.list = list(max.broad.roi = col.broi)



## Setting limit on max score to make colors easier to interpret
max = 10
go.mat[go.mat > max] = max
go.mat[go.mat < (-1*max)] = (-1*max)


## Plot heatmap
pheatmap(go.mat, fontsize_row = 6, fontsize_col = 6, #scale = "row",
         annotation_row = anno.df, annotation_colors = anno.list,
         filename = "figures_go/gprofiler_gomat.pdf", width = 15, height = 10, show_colnames = T,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         clustering_method = "ward.D", ##complete, ward.D, and ward.D2 are usually the best ones
         border_color = NA)