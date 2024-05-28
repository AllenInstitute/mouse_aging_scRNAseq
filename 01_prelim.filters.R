library(scrattch.hicat)
library(dplyr)

## Description/Goal:
## Filter data at the single-cell level using loose filter criteria

#################################################
## Load data files
#################################################
load("samp.dat_bothages_qcscore_20220705.rda") ## loading in cell metadata file

## Filter by very loose criteria
qc.cut = 50 ## qc score cutoff
gc.cut = 1000 ## gene detection cutoff
doub.cut = 0.3 ## doublet score cutoff

samp.dat.filtered = samp.dat[samp.dat$gene.counts.0 > gc.cut &
                               samp.dat$qc.score > qc.cut &
                               samp.dat$doublet_score < 0.3,]

save(samp.dat.filtered, file = "samp.dat.filtered_prelim.20220705.rda") ## save new filtered meta-data object for later