##########
# Project: Geographic potential of the world's largest hornet, Vespa mandarinia 
#          Smith (Hymenoptera: Vespidae), worldwide and particularly in North America
# Authors: 
# Claudia Nunez-Penichet, Luis Osorio-Olvera, Victor H. Gonzalez, Marlon E. Cobos, 
# Laura Jimenez, Devon A. DeRaad, Abdelghafar Alkishe, Rusby G. Contreras-Diaz, 
# Angela Nava-Bolanos, Kaera Utsumi, Uzma Ashraf, Adeola Adeboje, A. Townsend 
# Peterson, Jorge Soberon
#
# Date: 08/09/2020
##########

# Script to perform MOP analyses

# functions
source("R_code_V_mandarinia/00_Functions.R")

# directory
setwd("Working_directory")

#----------------------Performing MOP analyses accordingly----------------------
# MOPs for distance thinned records raw variables
## Reading variables to be used for MOPs
dtrv <- "Model_calibration/DT_RV"
mpath <- paste0(dtrv, "/all_results/03_mops_by_dataset/mops_table.csv")
mops <- read.csv(mpath)

## Identifying sets for MOPs
gsets <- id_mop_set(back_dir = paste0(dtrv, "/background"), mop_table = mops)

## MOPs
mopdir <- paste0(dtrv, "/MOP_reults")
dir.create(mopdir)
gsets <- unique(gsets)

mps <- lapply(1:length(gsets), function(y) {
  ## layers
  mvars <- read.csv(paste0(dtrv, "/Background/", gsets[y], ".csv"))
  gvars <- stack(list.files(paste0(dtrv, "/G_variables/", gsets[y], "/projection"), 
                            pattern = "asc$", full.names = T))
  
  ## mop and writing
  mop <- kuenm_mop(M.variables = mvars, G.stack = gvars, percent = 5, 
                   comp.each = 2000, parallel = F) # if have access to good computer use parallel = T
  writeRaster(mop, filename = paste0(mopdir, "/MOP_", gsets[y], "_5%.tif"), 
              format = "GTiff")
})


# MOPs for distance thinned records pcs
## Reading variables to be used for MOPs
dtpc <- "Model_calibration/DT_PC"
mpath <- paste0(dtpc, "/all_results/03_mops_by_dataset/mops_table.csv")
mops <- read.csv(mpath)

## Identifying sets for MOPs
gsets <- id_mop_set(back_dir = paste0(dtpc, "/background"), mop_table = mops)

## MOPs
mopdir <- paste0(dtpc, "/MOP_reults")
dir.create(mopdir)
gsets <- unique(gsets)

mps <- lapply(1:length(gsets), function(y) {
  ## layers
  mvars <- read.csv(paste0(dtpc, "/Background/", gsets[y], ".csv"))
  gvars <- stack(list.files(paste0(dtpc, "/G_variables/", gsets[y], "/projection"), 
                            pattern = "asc$", full.names = T))
  
  ## mop and writing
  mop <- kuenm_mop(M.variables = mvars, G.stack = gvars, percent = 5, 
                   comp.each = 2000, parallel = F) # if have access to good computer use parallel = T
  writeRaster(mop, filename = paste0(mopdir, "/MOP_", gsets[y], "_5%.tif"), 
              format = "GTiff")
})


# MOPs for country thinned records raw variables
## Reading variables to be used for MOPs
ctrv <- "Model_calibration/CT_RV"
mpath <- paste0(ctrv, "/all_results/03_mops_by_dataset/mops_table.csv")
mops <- read.csv(mpath)

## Identifying sets for MOPs
gsets <- id_mop_set(back_dir = paste0(ctrv, "/background"), mop_table = mops)

## MOPs
mopdir <- paste0(ctrv, "/MOP_reults")
dir.create(mopdir)
gsets <- unique(gsets)

mps <- lapply(1:length(gsets), function(y) {
  ## layers
  mvars <- read.csv(paste0(ctrv, "/Background/", gsets[y], ".csv"))
  gvars <- stack(list.files(paste0(ctrv, "/G_variables/", gsets[y], "/projection"), 
                            pattern = "asc$", full.names = T))
  
  ## mop and writing
  mop <- kuenm_mop(M.variables = mvars, G.stack = gvars, percent = 5, 
                   comp.each = 2000, parallel = F) # if have access to good computer use parallel = T
  writeRaster(mop, filename = paste0(mopdir, "/MOP_", gsets[y], "_5%.tif"), 
              format = "GTiff")
})


# MOPs for country thinned records pcs
## Reading variables to be used for MOPs
ctpc <- "Model_calibration/CT_PC"
mpath <- paste0(ctpc, "/all_results/03_mops_by_dataset/mops_table.csv")
mops <- read.csv(mpath)

## Identifying sets for MOPs
gsets <- id_mop_set(back_dir = paste0(ctpc, "/background"), mop_table = mops)

## MOPs
mopdir <- paste0(ctpc, "/MOP_reults")
dir.create(mopdir)
gsets <- unique(gsets)

mps <- lapply(1:length(gsets), function(y) {
  ## layers
  mvars <- read.csv(paste0(ctpc, "/Background/", gsets[y], ".csv"))
  gvars <- stack(list.files(paste0(ctpc, "/G_variables/", gsets[y], "/projection"), 
                            pattern = "asc$", full.names = T))
  
  ## mop and writing
  mop <- kuenm_mop(M.variables = mvars, G.stack = gvars, percent = 5, 
                   comp.each = 2000, parallel = F) # if have access to good computer use parallel = T
  writeRaster(mop, filename = paste0(mopdir, "/MOP_", gsets[y], "_5%.tif"), 
              format = "GTiff")
})
#-------------------------------------------------------------------------------
