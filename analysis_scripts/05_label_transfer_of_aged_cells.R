## Description/Goal:
## Performing label tansfer for all aged cells from ABC-WMB taxonomy (AIT21.0)

## Note:
## This script was run on an HPC job scheduler, and prior to the publishing of scrattch.mapping package
## We recomment installing github.com/alleninstitute/scrattch.mapping package to perform
## label transfer, but are including original scripts here for the sake of transparency

####################################################################
## These arguments were used for indexing on an HPC job scheduler ##
args = commandArgs(TRUE)
array.index = args[1]
array.index = as.integer(array.index)
####################################################################


## These functions have since been implemented in github.com/alleninstitute/scrattch-mapping
## but are also included in github repo for this project
source("functions/scrattch.mapping/HANN.R")
source("functions/scrattch.mapping/HANN_build.R")
source("functions/scrattch.mapping/HANN_utils.R")
source("functions/scrattch.mapping/HANN_prepareTaxonomy.R")


#################################################
## Load data files
#################################################
library(scrattch.bigcat)
library(openblasctl) 
openblas_set_num_threads(1)

load("big.dat_Aging.16ROI.V3.20231208.rda")
load("samp.dat_bothages_qcscore_20231208.rda")

## Select cells to keep
keep.cells = samp.dat$sample_id[samp.dat$age_cat == "aged"] ## Only select aged cells


#################################################
## Run mapping
#################################################
library(dplyr)
library(scrattch.bigcat)
library(Matrix)
library(data.table)
library(arrow)

## Breaking the data into parts for parallel processing on HPC
chunksize = 50000
start = seq(1,length(keep.cells), chunksize)
end = c(start[-1] - 1, length(keep.cells))
print(c(array.index, start[array.index], end[array.index]))


## Perform label transfer
select.cells = keep.cells[c(start[array.index]:end[array.index])]

print("Get norm.dat...")
qdat = get_cols(big.dat, select.cells)

print("Start mapping...")
mapped = run_mapping_on_taxonomy(
  qdat,                       # qdat : log normalized count matrix (gene x cell)
  Taxonomy='AIT21.0_mouse', # taxonomy
  TaxHome='/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/',
  prefix='10X_cells_v3',      # prefix for your data platform (SS, 10X_cells_v2, 10X_cells_v3, 10X_nuclei_v3,MFISH)
  prebuild=FALSE,             # if the taxonomy is available for your platform, use it.
  newbuild=FALSE,              #
  mapping.method='flat', # 'flat','hierarchy' (AIT12.0_mouse: only flat)
  nlevel=2,
  mc.cores=10,                # lower this to 5 if the run fails due to the memory
  iter=100,
  blocksize=5000)

fn.tmp = paste0("mapping_aged_fmap_AIT21.0_allcells/map.results_part", array.index, "of", length(start), ".rds")
print(paste0("Saving as ", fn.tmp, "..."))
saveRDS(mapped, file = fn.tmp)


#################################################
## Assemble mapping results
#################################################
fmap.files = list.files("./mapping_aged_fmap_AIT21.0/", pattern = "*.rds", full.names = T)
fmap.result = c()

for(i in fmap.files){
  
  tmp.fmap = readRDS(i)
  fmap.result = rbind(fmap.result, tmp.fmap$best.map.df)
  
}

head(fmap.result)
names(fmap.result)[-1] = paste0(names(fmap.result)[-1], "_AIT21.0")
ref.cl.df = tmp.fmap$cl.df

saveRDS(fmap.result, file = "mapping_aged_fmap_AIT21.0/map.result_bestmap.allcells.20230915.rds") ## Mapping results
save(ref.cl.df, file = "mapping_aged_fmap_AIT21.0/ref.cl.df.rda") ## Saves the annotation table for the reference that was used (ABC-WMB taxonomy)

