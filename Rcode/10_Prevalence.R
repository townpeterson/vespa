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

# Script for computing the prevalence of the species
# for all thinning and variable type scenarios:

# "50 km spatial thinning - PC",
# "50 km spatial thinning - raw variables",
# "Country-density thinning - PC variables",
# "Country-density thinning - raw variables"

# Packages
library(raster)
library(bam)
library(matrixStats)
library(dplyr)

# functions
source("R_code_V_mandarinia/00_Functions.R")

#-----------------Preparing data for prevalence calculations--------------------
# Folder path where some data are
dtrv <- "Model_calibration/DT_RV"

# Read presence data
native <- read.csv(list.files(dtrv, pattern = "all.csv$")[1])

# Read invasion points
invasion <- read.csv("coord_invasion.csv")

# Regex to look for extrapolative models projected in all the world
regxp <- "median_allSets_E.tif$"

# Paths of extrapolative models
models_p <- list.files(pattern = regxp,full.names = T,recursive = T)
models <- stack(models_p)

# Suitability values for native records 
nat_suits <- extract(models,native[,2:3])

# Suitability values for invasion records
inv_suits <- extract(models,invasion[,1:2])
inv_suits <- na.omit(inv_suits)
nat_suits <- na.omit(nat_suits)

# Suitability values for native and invasion records
all_suits <- rbind(nat_suits,inv_suits)

# Regex to look for extrapolative models projected in North America
regxpam <- "median_allSets_E_NA.tif$"

# Paths of extrapolative models in North America
models_amp <- list.files(pattern = regxpam,
                         full.names = T,
                         recursive = T)

# Read models
models_am <- stack(models_amp)
names(models_am) <- models_amp

# Five percent thresholds
th_vec <- c(0.05)
nlayers <- nlayers(models_am)

# Threshold values at five percent for each thinning scenario
rst_vals <- lapply(th_vec, function(x){
  mods_ras <- lapply(1:nlayers, function(y){
    r <- thfunc(median_mod = models_am[[y]],
                suits = all_suits[,y],percent = x)
  })
  return(mods_ras)
})

# Threshold values at five percent for each thinning scenario
suits_thresholds <- unlist(rst_vals[1])

# Binarize rasters 
modelos_bin <- models_am >= suits_thresholds
plot(modelos_bin)
#-------------------------------------------------------------------------------


#---------------------------Calculating prevalences-----------------------------
# Model names 
model_name <- c("50 km spatial thinning - PC",
                "50 km spatial thinning - raw variables",
                "Country-density thinning - PC variables",
                "Country-density thinning - raw variables")

# Look for pixels with values
nonasids <- which(!is.na(modelos_bin[[1]][]))

# Compute the area in km2 for Noth America
totalareaR <- area(modelos_bin[[1]],na.rm=T)

# Total area of North America
total_area <- sum(totalareaR[],na.rm = T)

# Chuck to compute the prevalance of models (without being trimmed by MOP)
rdf_nomop <- 1:length(model_name) %>% purrr::map_df(function(x){
  area2 <- totalareaR * modelos_bin[[x]]
  prev <- sum(area2[nonasids],na.rm = T)
  
  data.frame(model_name =model_name[x],
             total_area_km2=  total_area,
             predicted_presence_km2= prev,
             prevalence=prev/total_area,
             row.names = F,MOP="NO")
})                

row.names(rdf_nomop) <- NULL

# Read mops in North America
rasters_path <- list.files("MOPs_NA", pattern = "NA.tif$",recursive = T)
mops_na <- grep("MOP|mop",rasters_path,value = T)
mops_path <- c("50km_PC/abdu_na/sta.adf",
               mops_na)

# Filter non-extrapolative areas 
mops_bin <- stack(mops_path) == 0

# Trim models by MOPS
mod_bin_mop <- modelos_bin * mops_bin

# No NA ids
nonasids <- which(!is.na(mod_bin_mop[[1]][]))

# Chuck to compute the prevalance of models (trimmed by MOP)
rdf_mop <- 1:length(model_name) %>% purrr::map_df(function(x){
  area2 <- totalareaR * mod_bin_mop[[x]]
  prev <- sum(area2[nonasids],na.rm = T)
 
  data.frame(model_name =model_name[x],
             total_area_km2=  total_area,
             predicted_presence_km2= prev,
             prevalence=prev/total_area,
             row.names = F,MOP="YES")
})   

row.names(rdf_mop) <- NULL

# Join prevalence tables
tabla_conteos <- rbind.data.frame(rdf_mop,rdf_nomop)

# Write prevalence results
write.csv(tabla_conteos,"prevalence_km2.csv",row.names = F)
#-------------------------------------------------------------------------------
