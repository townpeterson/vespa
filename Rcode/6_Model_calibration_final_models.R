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

# Model calibration, Final model production 


# The next part helps in performing Model calibration, Final model production, 
# and model projections.
# 1) We will prepare the data the way it needs to be prepared (SWD format). This format
# does not use raster layers for calibration, instead we will use the values of such layers 
# extracted to csv files.
# 2) Model calibration will be done using the previously prepared data.
# 3) Final models will be constructed and projected to the World using the
# parameters selected during model calibration.


# install packages if needed
install.packages("remotes")

## kuenm needs Rtools or other compilation tools
## for more instructions see https://github.com/marlonecobos/kuenm#installing-the-package
remotes::install_github("marlonecobos/kuenm")

# load packages
library(raster)
library(kuenm)

# working directory to load initial data
setwd("working directory")

#-------------------------Preparing data for calibration------------------------
# number of model calibration exercises and variables
n <- 1:5

varsm <- stack(list.files("Model_calibration/Raw_variables_M", pattern = ".asc$",
                          full.names = T))

pcs <- stack(list.files("Model_calibration/PCs_M", pattern = ".asc$",
                        full.names = T))

# preparing data in loop for raw variables
## preparing data distance thinned records raw variables
dtrv <- "Model_calibration/DT_RV"
dir.create(dtrv)

alldt <- list.files("Model_calibration/Records_50km_thin/", pattern = "all.csv$")
traindt <- list.files("Model_calibration/Records_50km_thin/", pattern = "train.csv$")
testdt <- list.files("Model_calibration/Records_50km_thin/", pattern = "test.csv$")

ndt <- paste0(dtrv, "/vman")

prep <- lapply(n, function(w) {
  prep_swd(vars = varsm, all = read.csv(allr[w]), train = read.csv(trainr[w]), 
           test = read.csv(testdt[w]), name_occ = paste0(ndt, w),
           back_out_dir = paste0(dtrv, "/background"))
})


## preparing data distance thinned records PCs
dtpc <- "Model_calibration/DT_PC"
dir.create(dtpc)

alldt <- list.files("Model_calibration/Records_50km_thin/", pattern = "all.csv$")
traindt <- list.files("Model_calibration/Records_50km_thin/", pattern = "train.csv$")
testdt <- list.files("Model_calibration/Records_50km_thin/", pattern = "test.csv$")

ndt <- paste0(dtpc, "/vman")

prep <- lapply(n, function(w) {
  prep_swd(vars = pcs, all = read.csv(allr[w]), train = read.csv(trainr[w]), 
           test = read.csv(testdt[w]), name_occ = paste0(ndt, w),
           back_out_dir = paste0(dtpc, "/background"))
})


## preparing data country thinned records raw variables
ctrv <- "Model_calibration/CT_RV"
dir.create(ctrv)

alldt <- list.files("Model_calibration/Records_country_thin/", pattern = "all.csv$")
traindt <- list.files("Model_calibration/Records_country_thin/", pattern = "train.csv$")
testdt <- list.files("Model_calibration/Records_country_thin/", pattern = "test.csv$")

ndt <-  paste0(ctrv, "/vman")

prep <- lapply(n, function(w) {
  prep_swd(vars = varsm, all = read.csv(allr[w]), train = read.csv(trainr[w]), 
           test = read.csv(testdt[w]), name_occ = paste0(ndt, w),
           back_out_dir = paste0(ctrv, "/background"))
})


## preparing data country thinned records PCs
ctpc <- "Model_calibration/CT_PC"
dir.create(ctpc)

alldt <- list.files("Model_calibration/Records_country_thin/", pattern = "all.csv$")
traindt <- list.files("Model_calibration/Records_country_thin/", pattern = "train.csv$")
testdt <- list.files("Model_calibration/Records_country_thin/", pattern = "test.csv$")

ndt <- paste0(ctpc, "/vman")

prep <- lapply(n, function(w) {
  prep_swd(vars = pcs, all = read.csv(allr[w]), train = read.csv(trainr[w]), 
           test = read.csv(testdt[w]), name_occ = paste0(ndt, w),
           back_out_dir = paste0(ctpc, "/background"))
})
#-------------------------------------------------------------------------------


#-----------------------Model calibration and projections-----------------------
# general arguments
mx <- "C:/Maxent/3.4.1" # where is maxent in your computer (complete path)
worldraw <- "Model_calibration/Raw_variables_world"
worldpcs <- "Model_calibration/PCs_world"

# calibration and projections
## calibration and projection with data distance thinned records raw variables
all <- list.files(dtrv, pattern = "all.csv$")
train <- list.files(dtrv, pattern = "train.csv$")
test <- list.files(dtrv, pattern = "test.csv$")
bc <- paste0(dtrv, "/batch_cal", n)
cmdr <- paste0(dtrv, "/candidates", n)
calres <- paste0(dtrv, "/batch_cal", n)
bf <- paste0(dtrv, "/batch_fin", n)
fmods <- paste0(dtrv, "/final_models", n)

mod <- lapply(n, function(w) {
  maxent_calib_project(all_occ = all[w], tra_occ = train[w], tes_occ = test[w], 
                       back_dir = paste0(dtrv, "/background"), mx_path = mx, 
                       batch_cal = bc[w], candidate_models_dir = cmdr[w],
                       cal_res_dir = calres[w], all_gvars = worldraw, 
                       batch_fin = bf[w], final_models_dir = fmods[w])
})


## calibration and projection with data distance thinned records PCs
all <- list.files(dtpc, pattern = "all.csv$")
train <- list.files(dtpc, pattern = "train.csv$")
test <- list.files(dtpc, pattern = "test.csv$")
bc <- paste0(dtpc, "/batch_cal", n)
cmdr <- paste0(dtpc, "/candidates", n)
calres <- paste0(dtpc, "/batch_cal", n)
bf <- paste0(dtpc, "/batch_fin", n)
fmods <- paste0(dtpc, "/final_models", n)

mod <- lapply(n, function(w) {
  maxent_calib_project(all_occ = all[w], tra_occ = train[w], tes_occ = test[w], 
                       back_dir = paste0(dtpc, "/background"), mx_path = mx, 
                       batch_cal = bc[w], candidate_models_dir = cmdr[w],
                       cal_res_dir = calres[w], all_gvars = worldraw, 
                       batch_fin = bf[w], final_models_dir = fmods[w])
})


## calibration and projection with data country thinned records raw variables
all <- list.files(ctrv, pattern = "all.csv$")
train <- list.files(ctrv, pattern = "train.csv$")
test <- list.files(ctrv, pattern = "test.csv$")
bc <- paste0(ctrv, "/batch_cal", n)
cmdr <- paste0(ctrv, "/candidates", n)
calres <- paste0(ctrv, "/batch_cal", n)
bf <- paste0(ctrv, "/batch_fin", n)
fmods <- paste0(ctrv, "/final_models", n)

mod <- lapply(n, function(w) {
  maxent_calib_project(all_occ = all[w], tra_occ = train[w], tes_occ = test[w], 
                       back_dir = paste0(ctrv, "/background"), mx_path = mx, 
                       batch_cal = bc[w], candidate_models_dir = cmdr[w],
                       cal_res_dir = calres[w], all_gvars = worldraw, 
                       batch_fin = bf[w], final_models_dir = fmods[w])
})


## calibration and projection with data country thinned records PCs
all <- list.files(ctpc, pattern = "all.csv$")
train <- list.files(ctpc, pattern = "train.csv$")
test <- list.files(ctpc, pattern = "test.csv$")
bc <- paste0(ctpc, "/batch_cal", n)
cmdr <- paste0(ctpc, "/candidates", n)
calres <- paste0(ctpc, "/batch_cal", n)
bf <- paste0(ctpc, "/batch_fin", n)
fmods <- paste0(ctpc, "/final_models", n)

mod <- lapply(n, function(w) {
  maxent_calib_project(all_occ = all[w], tra_occ = train[w], tes_occ = test[w], 
                       back_dir = paste0(ctpc, "/background"), mx_path = mx, 
                       batch_cal = bc[w], candidate_models_dir = cmdr[w],
                       cal_res_dir = calres[w], all_gvars = worldraw, 
                       batch_fin = bf[w], final_models_dir = fmods[w])
})
#-------------------------------------------------------------------------------
