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

# Script for computing the binary maps of the selected MaxEnt models and
# perform simulations of the potential dynamics of dispersion for V. mandarinia

# Install the bam package 
# remotes::install_github("luismurao/bam")

# Loading packages
library(raster)
library(bam)
library(matrixStats)
library(animation)
library(furrr)
library(magrittr)

# functions
source("R_code_V_mandarinia/00_Functions.R")

# working directory
setwd("Working_directory")


#------------------------Preparing data for simulations-------------------------
# Folder path where some data are
dtrv <- "Model_calibration/DT_RV"

# Read presence data
native <- read.csv(list.files(dtrv, pattern = "all.csv$")[1])
invasion <- read.csv("coords_invasion_nortam.csv")

# Read SDM outputs
regxp <- "median_allSets_E.tif$"
models_p <- list.files(pattern = regxp,full.names = T,recursive = T)
models <- stack(models_p)

# Extract suitability values at presence points
nat_suits <- extract(models,native[,2:3])
inv_suits <- extract(models,invasion[,1:2])

# Remove NA values
inv_suits <- na.omit(inv_suits)
nat_suits <- na.omit(nat_suits)

# suitability values of invasion and native points
all_suits <- rbind(nat_suits,inv_suits)

# Path to models with free extrapolation projected in North America 
regxpam <- "median_allSets_E_NA.tif$"
models_amp <- list.files("Models_NA", pattern = regxpam,
                         full.names = T,
                         recursive = T)

# Create a raster stack from SDMs outputs
models_am <- stack(models_amp)
names(models_am) <- models_amp
nlayers <- nlayers(models_am)


# Define vector of threshold values
th_vec <- c(0.03,0.10)

# Apply threshold function for each threshold value and model
rst_vals <- lapply(th_vec, function(x){
  mods_ras <- lapply(1:nlayers, function(y){
    r <- thfunc(median_mod = models_am[[y]],
                suits = all_suits[,y],percent = x)
  })
  
  print(x)
  return(mods_ras)
})

# Suitability threshold values for 3%
suits_minimas <- unlist(rst_vals[1])

# Suitability threshold values for 10%
suits_maximas <- unlist(rst_vals[2])

# All suitability values in a matrix of 
# 2 (thresholds for 3% and 10%) X 4 (number of models representing the
#   spatial thinning scenario)
smin_max <-rbind(suits_minimas,suits_maximas)

# Binarize models using thresholds in the smin_max matrix 
modsbin <- lapply(1:ncol(smin_max), function(x){
  # A sequence of ten suitability values between 3% and 10% threshold
  umbrales <- seq(smin_max[1,x],
                  smin_max[2,x],
                  length.out = 10)
  
  # Binarize models
  r1 <- lapply(seq_along(umbrales), function(y){
    r0 <-models_am[[x]] >= umbrales[y]
  })
  names(r1) <- paste0("th_",round(umbrales,2))
  r1 <- stack(r1)
  return(r1)
})
#-------------------------------------------------------------------------------


#-----------Reducing suitable areas considering risks of extrapolation----------
# In this section we use the MOPs to trim the models
# corresponding to each spatial thinning scenario

# Read MOPs
rasters_path <- list.files("MOPs_NA", pattern = "NA.tif$",recursive = T)
mops_na <- grep("MOP|mop", rasters_path,value = T)
mops_path <- c("50km_PC/abdu_na/sta.adf",
               mops_na)

# Locating areas of extrapolation using the mops 
# This code generates binary MOPs where a value of 1
# means no extrapolation and 0 means extrapolation
mops_bin <- stack(mops_path) == 0

# Combine binary maps with MOPs for each scenario
# PC & 50 km 
maps_50km_PC <- modsbin[[1]] * mops_bin[[1]]

# Raw & 50 km 
maps_50km_rw <- modsbin[[2]] * mops_bin[[2]]

# PC & country 
maps_ct_PC <- modsbin[[3]] * mops_bin[[3]]

# Raw & country 
maps_ct_rw <- modsbin[[4]] * mops_bin[[4]]

# Plot the models for 
# 3% threshold and tenth percentile for each scenario
# PC & 50 km
plot(maps_50km_PC[[c(1,10)]])

# Raw & 50 km 
plot(maps_50km_rw [[c(1,10)]])

# PC & country 
plot(maps_ct_PC[[c(1,10)]])

# Raw & country 
plot(maps_ct_rw[[c(1,10)]])
#-------------------------------------------------------------------------------


#------------------------------Dispersal simulations----------------------------
# Connectivity matrices
# first, we convert the model into a diagonal sparse matrix
sparse_temp <-  bam::model2sparse(maps_50km_PC[[1]])

# and, we convert the occurrence points into a sparse matrix
occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp,
                                occs = invasion[,1:2])

# second, we create adjacency matrices for different distance values 
ad_L <- lapply(c(1:8,10,12), function(x){
  adj_matrices <- bam::adj_mat(sparse_temp,ngbs=x)
  print(x)
  return(adj_matrices)
})

# Names of the matrices 
names(adL) <- c(1:8,10,12)

# Third, determine the number of steps to be run in the simulation process
pasos <- 200

# Finally, apply the simulation function to each scenario
# PC & 50 km 
r_50km_PC <- sim_disperal(adL = ad_L, bin_thmodels = maps_50km_PC,
                          base_name = paste0("simul_res_50km_PC_",pasos))

# Raw & 50 km 
r_50km_rw <- sim_disperal(adL = ad_L, bin_thmodels = maps_50km_rw,
                          base_name = paste0("simul_res_50km_rw_",pasos))

# PC & country 
r_ct_PC <- sim_disperal(adL = ad_L, bin_thmodels = maps_ct_PC,
                          base_name = paste0("simul_res_ct_PC_",pasos))

# Raw & country
r_ct_rw <- sim_disperal(adL = ad_L, bin_thmodels = maps_ct_rw,
                          base_name = paste0("simul_res_ct_rw_",pasos))


# Calculate the sum of all the resulting maps for the different distance values
r_50km_PC_sum <- calc(stack(r_50km_PC),sum)
plot(r_50km_PC_sum)

r_50km_rw_sum <- calc(stack(r_50km_rw),sum)
plot(r_50km_rw_sum)

r_ct_PC_sum <- calc(stack(r_ct_PC),sum)
plot(r_ct_PC_sum)

r_ct_rw_sum <- calc(stack(r_ct_rw),sum)
plot(r_ct_rw_sum)

# save results
writeRaster(r_50km_PC_sum,paste0("r_50km_PC_sum_",pasos,".tif"),format="GTiff", 
            overwrite=TRUE)

writeRaster(r_50km_rw_sum,paste0("r_50km_rw_sum_",pasos,".tif"),format="GTiff", 
            overwrite=TRUE)

writeRaster(r_ct_PC_sum,paste0("r_ct_PC_sum_",pasos,".tif"),format="GTiff", 
            overwrite=TRUE)

writeRaster(r_ct_rw_sum,paste0("r_ct_rw_sum_",pasos,".tif"),format="GTiff", 
            overwrite=TRUE)
#-------------------------------------------------------------------------------

