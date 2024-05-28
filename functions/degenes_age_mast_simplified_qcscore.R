get_age_degenes_mast = function(anno.df, norm.dat, id, split.by, maxcells = 500, freq_expressed = 0.1, exclude = F, save.as.tmp = F, select.genes = NULL, select.cells = NULL){
	
	set.seed(2010)

	require(MAST)
	require(dplyr)
	require(scrattch.hicat)
	require(data.table)

	## Rename objects
	anno.df.select = anno.df
	row.names(anno.df.select) = anno.df.select$sample_id
	group.tmp = id
	dat.select = norm.dat
	print(paste0("Starting ", group.tmp, "..."))


    ## Exclude or include
    if(exclude){
    	anno.tmp = anno.df.select[anno.df.select[,split.by] != group.tmp,]
    } else {
    	anno.tmp = anno.df.select[anno.df.select[,split.by] == group.tmp,]

    }

    ## Subsample to the same size per age group, up to maxcells
    if(is.null(select.cells)){
      tmp.age = anno.tmp$sample_id[anno.tmp$age_cat == "aged"]
      tmp.adult = anno.tmp$sample_id[anno.tmp$age_cat == "adult"]
      tmp.max.age = min(length(tmp.age), maxcells)
      tmp.max.adult = min(length(tmp.adult), maxcells)
      select.cells = c(sample(tmp.age, tmp.max.age), sample(tmp.adult, tmp.max.adult))
      anno.tmp = anno.tmp[select.cells,]
      anno.tmp$sex = as.character(anno.tmp$sex) ## correcting weird error where some entries were boolean rather than characters

      ## Escape for error
      tmp = c(table(anno.tmp[,c("sex", "age_cat")]))
      if(sum(tmp == 0) > 0){
        stop("Error: Missing data for one of the sexes and/or ages.")
        # break
      }
    }


    ## Subset data
    dat.tmp = dat.select[,select.cells]
    anno.tmp = anno.tmp[select.cells,]
    

    ## Make MAST object
    row.names(anno.tmp) = anno.tmp$sample_id
    sca = FromMatrix(as.matrix(dat.tmp), cData = anno.tmp)

  # ## Thresholding ##SKIPPING FOR NOW
  # thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 30, min_per_bin = 30)
  # par(mfrow=c(5,4))
  # plot(thres)

    ## Freq cutoff
    if(is.null(select.genes)){
      select.genes = freq(sca) > freq_expressed

    }
    print(paste0("Ncells: ", length(select.cells)))
    print(paste0("Ngenes: ", sum(select.genes)))
    sca = sca[select.genes,]
    
    ## Model
    print("Run model")
    colData(sca)$age_cat<-factor(colData(sca)$age_cat)
    colData(sca)$sex<-factor(colData(sca)$sex)
    colData(sca)$roi<-factor(colData(sca)$roi)
    colData(sca)$broad_roi<-factor(colData(sca)$broad_roi)
    colData(sca)$facs_population_plan<-factor(colData(sca)$facs_population_plan)
    colData(sca)$full_genotype<-factor(colData(sca)$full_genotype)
    colData(sca)$log.gene.counts.0 = scale(log(colData(sca)$gene.counts.0))
    colData(sca)$z.qc.score = scale(colData(sca)$qc.score)
    
    zlm <- zlm(~age_cat + sex + log.gene.counts.0 + z.qc.score, sca)


    ## ----> LTR
    summaryAge<- summary(zlm, doLRT='age_cataged')
    summaryDt.age <- summaryAge$datatable
    fcHurdle.age <- merge(summaryDt.age[contrast=='age_cataged' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt.age[contrast=='age_cataged' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle.age[,padjust:=p.adjust(`Pr(>Chisq)`, 'bonferroni')]
    fcHurdle.age[split.by] = group.tmp
    
    ## ----> Save results
    result = list(select.cells = select.cells, result = fcHurdle.age)
    if(save.as.tmp){
      save(result, file = paste0("tmp.degene.result.", make.names(id), ".rda"))
    }
    return(result)

}





get_age_degenes_mast_big = function(anno.df, big.dat, id, split.by, maxcells = 500, freq_expressed = 0.1, exclude = F, save.as.tmp = F, select.genes = NULL, select.cells = NULL){
	
	set.seed(2010)

	require(MAST)
	require(dplyr)
	require(scrattch.hicat)
	require(data.table)

	## Rename objects
	anno.df.select = anno.df
	row.names(anno.df.select) = anno.df.select$sample_id
	group.tmp = id
	print(paste0("Starting ", group.tmp, "..."))


    ## Exclude or include
    if(exclude){
    	anno.tmp = anno.df.select[anno.df.select[,split.by] != group.tmp,]
    } else {
    	anno.tmp = anno.df.select[anno.df.select[,split.by] == group.tmp,]

    }

    ## Subsample to the same size per age group, up to maxcells (if no preselected cells)
    if(is.null(select.cells)){
      tmp.age = anno.tmp$sample_id[anno.tmp$age_cat == "aged"]
      tmp.adult = anno.tmp$sample_id[anno.tmp$age_cat == "adult"]
      tmp.max.age = min(length(tmp.age), maxcells)
      tmp.max.adult = min(length(tmp.adult), maxcells)
      select.cells = c(sample(tmp.age, tmp.max.age), sample(tmp.adult, tmp.max.adult))
      anno.tmp = anno.tmp[select.cells,]
      anno.tmp$sex = as.character(anno.tmp$sex) ## correcting weird error where some entries were boolean rather than characters

      ## Escape for error 
      tmp = c(table(anno.tmp[,c("sex", "age_cat")]))
      if(sum(tmp == 0) > 0){
        stop("Error: Missing data for one of the sexes and/or ages.")
        # break
      }
    }


    ## Subsample data
    dat.tmp = get_logNormal(big.dat, select.cells)

    ## Make MAST object
    row.names(anno.tmp) = anno.tmp$sample_id
    sca = FromMatrix(as.matrix(dat.tmp), cData = anno.tmp)

  # ## Thresholding ##SKIPPING FOR NOW
  # thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 30, min_per_bin = 30)
  # par(mfrow=c(5,4))
  # plot(thres)

    ## Freq cutoff
    if(is.null(select.genes)){
      select.genes = freq(sca) > freq_expressed
    }
    print(paste0("Ncells: ", length(select.cells)))
    print(paste0("Ngenes: ", sum(select.genes)))
    sca = sca[select.genes,]
    
    ## Model
    print("Run model")
    colData(sca)$age_cat<-factor(colData(sca)$age_cat)
    colData(sca)$sex<-factor(colData(sca)$sex)
    colData(sca)$roi<-factor(colData(sca)$roi)
    colData(sca)$broad_roi<-factor(colData(sca)$broad_roi)
    colData(sca)$facs_population_plan<-factor(colData(sca)$facs_population_plan)
    colData(sca)$full_genotype<-factor(colData(sca)$full_genotype)
    colData(sca)$log.gene.counts.0 = scale(log(colData(sca)$gene.counts.0))
    colData(sca)$z.qc.score = scale(colData(sca)$qc.score)
    
    zlm <- zlm(~age_cat + sex + log.gene.counts.0 + z.qc.score, sca)
    


    ## ----> LTR
    summaryAge<- summary(zlm, doLRT='age_cataged')
    summaryDt.age <- summaryAge$datatable
    fcHurdle.age <- merge(summaryDt.age[contrast=='age_cataged' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt.age[contrast=='age_cataged' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle.age[,padjust:=p.adjust(`Pr(>Chisq)`, 'bonferroni')]
    fcHurdle.age[split.by] = group.tmp
    
    ## ----> Save results
    result = list(select.cells = select.cells, result = fcHurdle.age)
    if(save.as.tmp){
    	save(result, file = paste0("tmp.degene.result.", make.names(id), ".rda"))
    }
    return(result)

}


# 
# 
# 
# 
# 
# get_sex_intergenes_mast_big = function(anno.df, big.dat, id, split.by, maxcells = 500, freq_expressed = 0.1, exclude = F, save.as.tmp = F, select.genes = NULL, select.cells = NULL){
#   
#   set.seed(2010)
#   
#   require(MAST)
#   require(dplyr)
#   require(scrattch.hicat)
#   require(data.table)
#   
#   ## Rename objects
#   anno.df.select = anno.df
#   row.names(anno.df.select) = anno.df.select$sample_id
#   group.tmp = id
#   print(paste0("Starting ", group.tmp, "..."))
#   
#   
#   ## Exclude or include
#   if(exclude){
#     anno.tmp = anno.df.select[anno.df.select[,split.by] != group.tmp,]
#   } else {
#     anno.tmp = anno.df.select[anno.df.select[,split.by] == group.tmp,]
#     
#   }
#   
#   ## Subsample to the same size per age group, up to maxcells (if no preselected cells)
#   if(is.null(select.cells)){
#     tmp.age = anno.tmp$sample_id[anno.tmp$age_cat == "aged"]
#     tmp.adult = anno.tmp$sample_id[anno.tmp$age_cat == "adult"]
#     tmp.max.age = min(length(tmp.age), maxcells)
#     tmp.max.adult = min(length(tmp.adult), maxcells)
#     select.cells = c(sample(tmp.age, tmp.max.age), sample(tmp.adult, tmp.max.adult))
#     anno.tmp = anno.tmp[select.cells,]
#     anno.tmp$sex = as.character(anno.tmp$sex) ## correcting weird error where some entries were boolean rather than characters
#     
#     ## Escape for error 
#     tmp = c(table(anno.tmp[,c("sex", "age_cat")]))
#     if(sum(tmp == 0) > 0){
#       stop("Error: Missing data for one of the sexes and/or ages.")
#       # break
#     }
#   }
#   
#   
#   ## Subsample data
#   dat.tmp = get_logNormal(big.dat, select.cells)
#   
#   ## Make MAST object
#   row.names(anno.tmp) = anno.tmp$sample_id
#   sca = FromMatrix(as.matrix(dat.tmp), cData = anno.tmp)
#   
#   # ## Thresholding ##SKIPPING FOR NOW
#   # thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 30, min_per_bin = 30)
#   # par(mfrow=c(5,4))
#   # plot(thres)
#   
#   ## Freq cutoff
#   if(is.null(select.genes)){
#     select.genes = freq(sca) > freq_expressed
#   }
#   print(paste0("Ncells: ", length(select.cells)))
#   print(paste0("Ngenes: ", sum(select.genes)))
#   sca = sca[select.genes,]
#   
#   ## Model
#   print("Run model")
#   colData(sca)$age_cat<-factor(colData(sca)$age_cat)
#   colData(sca)$sex<-factor(colData(sca)$sex)
#   colData(sca)$roi<-factor(colData(sca)$roi)
#   colData(sca)$facs_population_plan<-factor(colData(sca)$facs_population_plan)
#   colData(sca)$full_genotype<-factor(colData(sca)$full_genotype)
#   
#   
#   
#   
#   
#   if(length(unique(anno.tmp$full_genotype)) == 1){
#     if(length(unique(anno.tmp$facs_population_plan)) == 1){
#       if(length(unique(anno.tmp$roi)) == 1){
#         zlm <- zlm(~age_cat + sex + age_cat:sex, sca)
#       } else{
#         zlm <- zlm(~age_cat + sex + roi + age_cat:sex, sca)
#       }
#     } else{
#       if(length(unique(anno.tmp$roi)) == 1){
#         zlm <- zlm(~age_cat + sex + facs_population_plan + age_cat:sex, sca)
#       } else{
#         zlm <- zlm(~age_cat + sex + roi + facs_population_plan + age_cat:sex, sca)
#       }
#     }
#   } else{
#     if(length(unique(anno.tmp$facs_population_plan)) == 1){
#       if(length(unique(anno.tmp$roi)) == 1){
#         zlm <- zlm(~age_cat + sex + full_genotype + age_cat:sex, sca)
#       } else{
#         zlm <- zlm(~age_cat + sex + roi + full_genotype + age_cat:sex, sca)
#       }
#     } else{
#       if(length(unique(anno.tmp$roi)) == 1){
#         zlm <- zlm(~age_cat + sex + facs_population_plan + full_genotype + age_cat:sex, sca)
#       } else{
#         zlm <- zlm(~age_cat + sex + roi + facs_population_plan + full_genotype + age_cat:sex, sca)
#       }
#     }
#   }
#   
#   
#   
#   
#   ## ----> LTR
#   summaryAge<- summary(zlm, doLRT='age_cataged:sexM')
#   summaryDt.age <- summaryAge$datatable
#   fcHurdle.age <- merge(summaryDt.age[contrast=='age_cataged:sexM' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
#                         summaryDt.age[contrast=='age_cataged:sexM' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
#   fcHurdle.age[,padjust:=p.adjust(`Pr(>Chisq)`, 'bonferroni')]
#   fcHurdle.age[split.by] = group.tmp
#   
#   ## ----> Save results
#   result = list(select.cells = select.cells, result = fcHurdle.age)
#   if(save.as.tmp){
#     save(result, file = paste0("tmp.degene.result.", make.names(id), ".rda"))
#   }
#   return(result)
#   
# }








get_age_degenes_mast_p = function(anno.df, norm.dat, groups, split.by, maxcells = 500, freq_expressed = 0.1, mc.cores = 10, exclude = F, save.as.tmp = F, select.genes = NULL, select.cells = NULL){

	if (mc.cores == 1) {
        
        de.list = list()

        for (n in 1:length(groups)) {

        	id = groups[n]
            de.list[id] = list(get_age_degenes_mast(anno.df, norm.dat, id = id, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = select.cells))

        }

    }else {

	    library(doMC)
	    registerDoMC(cores = mc.cores)
	    de.list = foreach::foreach(i = groups, .combine='c') %dopar% {
	                 list(get_age_degenes_mast(anno.df, norm.dat, id = i, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = select.cells))
	              }
	    names(de.list) = groups


    }

    return(de.list)

}







get_age_degenes_mast_big_p = function(anno.df, big.dat, groups, split.by, maxcells = 500, freq_expressed = 0.1, mc.cores = 10, exclude = F, save.as.tmp = F, select.genes = NULL, select.cells = NULL){

	if (mc.cores == 1) {
        
        de.list = list()

        for (n in 1:length(groups)) {

        	id = groups[n]
            de.list[id] = list(get_age_degenes_mast_big(anno.df, big.dat, id = id, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = select.cells))

        }

    }else {

	    library(doMC)
	    registerDoMC(cores = mc.cores)
	    de.list = foreach::foreach(i = groups, .combine='c') %dopar% {
	                 list(get_age_degenes_mast_big(anno.df, big.dat, id = i, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = select.cells))
	              }
	    names(de.list) = groups


    }

    return(de.list)

}



# 
# 
# get_sex_intergenes_mast_big_p = function(anno.df, big.dat, groups, split.by, maxcells = 500, freq_expressed = 0.1, mc.cores = 10, exclude = F, save.as.tmp = F, select.genes = NULL, select.cells = NULL){
#   
#   if (mc.cores == 1) {
#     
#     de.list = list()
#     
#     for (n in 1:length(groups)) {
#       
#       id = groups[n]
#       de.list[id] = list(get_sex_intergenes_mast_big(anno.df, big.dat, id = id, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = select.cells))
#       
#     }
#     
#   }else {
#     
#     library(doMC)
#     registerDoMC(cores = mc.cores)
#     de.list = foreach::foreach(i = groups, .combine='c') %dopar% {
#       list(get_sex_intergenes_mast_big(anno.df, big.dat, id = i, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = select.cells))
#     }
#     names(de.list) = groups
#     
#     
#   }
#   
#   return(de.list)
#   
# }
# 
# 







get_age_degenes_perm = function(anno.df, norm.dat, id, split.by, nperm = 50, maxcells = 500, freq_expressed = 0.1, mc.cores = 10, exclude = F, save.as.tmp = F, select.genes = NULL){
  
  select.cells.age = anno.df$sample_id[anno.df$age_cat == "aged"]
  select.cells.adult = anno.df$sample_id[anno.df$age_cat == "adult"]
  
  set.seed(2010)
  
  perm.mat.age = replicate(nperm, sample(select.cells.age, maxcells, replace = F))
  perm.mat.adult = replicate(nperm, sample(select.cells.adult, maxcells, replace = F))
  perm.mat = rbind(perm.mat.age, perm.mat.adult)
  
  if (mc.cores == 1) {
    
    de.list = list()
    
    for (n in 1:nperm) {
      
      pname = paste0("p",n)
      de.list[pname] = list(get_age_degenes_mast(anno.df, norm.dat, id = id, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = perm.mat[,n]))
      
    }
    
  }else {
    
    library(doMC)
    registerDoMC(cores = mc.cores)
    de.list = foreach::foreach(i = 1:nperm, .combine='c') %dopar% {
      list(get_age_degenes_mast(anno.df, norm.dat, id = id, split.by = split.by, maxcells = maxcells, freq_expressed = freq_expressed, exclude = exclude, save.as.tmp = save.as.tmp, select.genes = select.genes, select.cells = perm.mat[,i]))
    }
    names(de.list) = paste0("p",1:nperm)
    
    
  }
  
  return(de.list)
  
}
