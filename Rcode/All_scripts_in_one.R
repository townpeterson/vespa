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


#-------------------Installing and loading packages-----------------------------
#installing if needed
#install.packages("spocc")
#install.packages("rgbif")
#install.packages("maps")
#install.packages("rgdal")
#install.packages("rgeos")
#install.packages("raster")
#install.packages("sp")
#install.packages("dplyr")
#install.packages("purrr")
#install.packages("matrixStats")
#install.packages("animation")
#install.packages("furrr")
#install.packages("magrittr")
#
## url https://github.com/marlonecobos/ellipsenm#installing-the-package
#remotes::install_github("marlonecobos/ellipsenm")
## url https://github.com/marlonecobos/kuenm#installing-the-package
#remotes::install_github("marlonecobos/kuenm")
## url https://github.com/luismurao/ntbox#complete-installation-guide
#remotes::install_github("luismurao/ntbox")
## url https://github.com/luismurao/bam#installation
# remotes::install_github("luismurao/bam")

# packages
library(spocc)
library(rgbif)
library(maps)
library(ellipsenm)
library(rgdal)
library(rgeos)
library(raster)
library(dplyr)
library(purrr)
library(kuenm)
library(ntbox)
library(bam)
library(matrixStats)
library(animation)
library(furrr)
library(magrittr)
#-------------------------------------------------------------------------------


#-------------------Bringing functions and setting directory--------------------
# functions
source("R_code_V_mandarinia/00_Functions.R")

# directory
setwd("Working_directory")
#-------------------------------------------------------------------------------



#-------------------Downloading species occurrences-----------------------------
# getting occurrences from GBIF
vm <- occ("Vespa mandarinia", limit = 10000)
occs <- vm$gbif$data$Vespa_mandarinia
head(occs)

## references
sp_search <- occ_search(taxonKey = vm$gbif$data$Vespa_mandarinia$taxonKey[1])
cit <- gbif_citation(sp_search)
sink("gbif_ref.txt")
sapply(cit, print)
sink()

## saving initial data with various columns from GBIF
colnames(occs)
columns <- c("name", "longitude", "latitude", "issues", "scientificName", 
             "coordinateUncertaintyInMeters", "year", "month", "day", 
             "countryCode", "locality", "elevation")

occs <- occs[, columns]

write.csv(occs, "occurrences_GBIF.csv", row.names = F)
#-------------------------------------------------------------------------------


#-------------------- ----Data cleaning-----------------------------------------
colnames(occs)

## only columns of interest
occ <- occs[, c("scientificName", "longitude", "latitude", "year")]
occ <- na.omit(occ)

## defining best time interval for variables
hist(occ$year, breaks = 20)

## excluding older records
occ <- occ[occ$year >= 1990, -4]

## changing names in first column
unique(occ$scientificName)
occ <- occ[occ$scientificName == "Vespa mandarinia Smith, 1852", ]
occ$scientificName <- "Vespa mandarinia"

## excluding 0, 0 coordinates 
occ <- occ[occ$longitude != 0 & occ$latitude != 0, ]

## excluding duplicates
occ <- occ[!duplicated(paste(occ$longitude, occ$latitude)), ]

map()
points(occ[, 2:3], col = "red", pch = 19)
axis(side = 1)

## excluding record in Europe
occ <- occ[occ$longitude > 30, ]

map()
points(occ[, 2:3], col = "red", pch = 19)

write.csv(occ, "V_mandarinia_clean.csv", row.names = F)
#-------------------------------------------------------------------------------


#-----------------------------Thinning------------------------------------------
# spatial distance thinning
occt <- thin_data(occ, "longitude", "latitude", thin_distance = 50, save = T, 
                  name = "V_mandarinia_50km")

map(xlim = range(occt$longitude), ylim = range(occt$latitude))
points(occt[, 2:3], col = "red", pch = 19)


# country density thinning
set.seed(111)

recs_by_country <- c("Taiwan" = 2, "Japan" = 6, "South Korea" = 2)

# reading shapefile (the world shape file used here was: ESRI UIA_World Countries Boundaries, 
#2020. https://www.arcgis.com/home/item.html?id=252471276c9941729543be8789e06e12)
wrld <- readOGR("WORLD","country")
plot(wrld)

# extracting country names by coordinates in our data base
# Convert data.frame to SpatialPointsDataFrame
vm_sp <- vm
names(vm_sp)
coordinates(vm_sp) <- ~longitude + latitude
crs(vm_sp) <- crs(wrld)
plot(vm_sp, add = T)

# Get the shapefile attributes 
# Extract values to points 
vm_atributtes <- over(vm_sp, wrld)

# Join the attributes 
vm_df <- data.frame(vm, vm_atributtes)
vm_df$CNTRY_NAME <- as.character(vm_df$CNTRY_NAME)
hist(table(vm_df$CNTRY_NAME))

# Split data by country 
cto_sampleL <- vm_df %>% split(.$CNTRY_NAME,drop = T)
x = 2

samp_countries <- 1:length(cto_sampleL) %>%
  purrr::map_df(function(x){
    if(cto_sampleL[[x]]$CNTRY_NAME[1] == "Japan"){
      df <- cto_sampleL[[x]] %>% sample_n(6)
      return(df)
      
    }
    if(cto_sampleL[[x]]$CNTRY_NAME[1] == "South Korea"){
      df <- cto_sampleL[[x]] %>% sample_n(2)
      return(df)
    }
    if(cto_sampleL[[x]]$CNTRY_NAME[1] == "Taiwan"){
      df <- cto_sampleL[[x]] %>% sample_n(2)
      return(df)
      
    }
    else
      df <- cto_sampleL[[x]]
    return(df)
  })

points(samp_countries[, c(2, 3)], col = )

write.csv(samp_countries,"V_mandarinia_country.csv", row.names = F)

#-------------------------------------------------------------------------------


#-----------------Splitting training and testing data---------------------------
# model calibration folder
dir.create("Model_calibration")
dir.create("Model_calibration/Records_50km_thin")
dir.create("Model_calibration/Records_country_thin")

# if require to use the data provided, use the following lines
#occt <- read.csv("V_mandarinia_50km.csv)
#samp_countries <- read.csv("V_mandarinia_country.csv)

# split distance based thinned data for 5 exercises of model calibration
n <- 1:5
splits <- lapply(n, function(x) {
  set.seed(x)
  occsp_50 <- split_data(occt, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("Model_calibration/Records_50km_thin/vman", x))
})

# split country density data for 5 exercises of model calibration
splits <- lapply(n, function(x) {
  set.seed(x)
  occsp_50 <- split_data(samp_countries, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.5, save = T, 
                         name = paste0("Model_calibration/Records_country_thin/vman", x))
})
#-------------------------------------------------------------------------------


#-------------------Downloading, unzipping variables----------------------------
## download url
past_var <- "https://datadryad.org/stash/downloads/file_stream/97895"

# Function to download 
download.file(past_var, destfile = "10m_mean_00s.zip", method = "auto", 
              quiet = FALSE, mode = "wb", cacheOK = TRUE)

# unzipping variables
dfol <- "10m_mean_00s"
dir.create(dfol)
unzip(zipfile = "10m_mean_00s.zip", exdir = dfol)
#-------------------------------------------------------------------------------


#----------------------------Preparing variables--------------------------------
# stacking of variables (excluding 8, 9, 18, 19)
mc10 <- stack(list.files("10m_mean_00s", pattern = "bio.*if$", full.names = T))
names(mc10) <- gsub("X10m_mean_00s_", "", names(mc10))
mc10 <- mc10[[c(1, 12:17, 2:9)]]

# Problem: we have Sea temperatures but the species is terrestrial
# We need to set sea data as NA (the world shape file used here was: 
# ESRI UIA_World Countries Boundaries, 2020. 
# https://www.arcgis.com/home/item.html?id=252471276c9941729543be8789e06e12)
wrld <- readOGR("WORLD", "country")
wrld$id <- "1"

# Dissolve the polygon
nwd <- gUnaryUnion(wrld, wrld@data$id)

# Getting the cells with values 
cells_nona <- cellFromPolygon(mc10[[1]], nwd )
cellsids <- 1: ncell(mc10[[1]])
cell_na <- cellsids[-cells_nona[[1]]]
mc10[cell_na] <- NA
plot(mc10[[1]])

# Saving variables to new directory
dir.create("merra_new")

lapply(names(mc10), function(x){
  writeRaster(mc10[[x]],paste0("merra_new/", x, ".asc"), overwrite = T)
})
#-------------------------------------------------------------------------------


#-------------------------Preparing calibration area----------------------------
# occurrence data
occ <- read.csv("V_mandarinia_50km.csv")

# creating M
## considering earth distortion
WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
occ_sp <- SpatialPointsDataFrame(coords = occ[, 2:3], data = occ,
                                 proj4string = WGS84)

## project the points using their centroids as reference
centroid <- gCentroid(occ_sp, byid = FALSE)
AEQD <- CRS(paste("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", centroid@coords[1],
                  " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", sep = ""))

occ_pr <- spTransform(occ_sp, AEQD)

## create a buffer based on a 500 km distance
buff_area <- gBuffer(occ_pr, width = 500000, quadsegs = 30)
buff_area <- disaggregate(buff_area)

## reproject
buff_area <- spTransform(buff_area, WGS84)

## make spatialpolygondataframe
df <- data.frame(species = rep("V_mandarinia", length(buff_area)))
buff_area <- SpatialPolygonsDataFrame(buff_area, data = df, match.ID = FALSE)

## write area as shapefile
dir.create("calibration_area")
writeOGR(buff_area, "calibration_area", "M", driver = "ESRI Shapefile")

## plot
lims <- extent(buff_area)
map(xlim = lims[1:2] + c(-1, 1), ylim = lims[3:4] + c(-1, 1))
points(occ[, 2:3], col = "red", pch = 19)
plot(buff_area, border = "purple", add = T, lwd = 2)
#-------------------------------------------------------------------------------


#---------------------Masking variables to calibration area---------------------
# renaming M 
M <- buff_area

# masking layers to M
varsm <- mask(crop(mc10, M), M)

# saving masked layers as ascii files
lapply(names(varsm), function(x){
  writeRaster(varsm[[x]], paste0("calibration_area/", x,".asc"), overwrite=T)
})
#-------------------------------------------------------------------------------


#-------------------Exploring variables in environmental space------------------
# exploring the data in environmental space
## only temperature variables
explore_espace(data = occt, species = "scientificName", longitude = "longitude",
               latitude = "latitude", raster_layers = varsm[[1:9]], save = T,
               name = "Temperature_variables.pdf")

## only precipitation variables
explore_espace(data = occt, species = "scientificName", longitude = "longitude",
               latitude = "latitude", raster_layers = varsm[[10:15]], save = T,
               name = "Precipitation_variables.pdf")

## highly correlated set
hcor <- c(1, 6, 8, 9, 10, 11, 12, 14, 15)
explore_espace(data = occt, species = "scientificName", longitude = "longitude",
               latitude = "latitude", raster_layers = varsm[[hcor]], 
               save = T, name = "High-cor_variables.pdf")

## non-correlated sets
ncor_sets <- list(c(3, 5, 7, 11, 12), c(3, 4, 5, 11, 12), c(2, 3, 10, 13), 
                  c(3, 4, 10, 13), c(3, 4, 12, 13), c(3, 7, 12, 13), 
                  c(3, 4, 13, 15), c(3, 7, 13, 15))

explore <- lapply(3:length(ncor_sets), function(x) {
  explore_espace(data = occt, species = "scientificName", longitude = "longitude",
                 latitude = "latitude", raster_layers = varsm[[ncor_sets[[x]]]], 
                 save = T, name = paste0("Non-cor_variables", x, ".pdf"))
})


# exploring variable correlation in one plot for all
jpeg("corrplot_merra.jpg", width = 120, height = 120, units = "mm", res = 600)
par(cex = 0.8)
vcor <- variable_correlation(varsm, save = T, name = "correlation_merra", 
                             corrplot = T, magnify_to = 3)
dev.off()

# variables selected were bio: 3, 5, 6, 7, 13, 14
dir.create("Model_calibration/Raw_variables_M")

nums <- c(3, 5, 6, 7, 13, 14)

file.copy(from = paste0("calibration_area/bio", nums,".asc"), 
          to = paste0("Model_calibration/Raw_variables_M", nums,".asc"))

dir.create("Model_calibration/Raw_variables_world")

file.copy(from = paste0("merra_new/bio", nums,".asc"), 
          to = paste0("Model_calibration/Raw_variables_world", nums,".asc"))
#-------------------------------------------------------------------------------


#------------------Principal component analysis and projections-----------------
# PCA and projections
dir.create("pcas")
dir.create("pcas/pca_referenceLayers")
dir.create("pcas/pca_proj")
s1 <- spca(layers_stack = varsm, layers_to_proj = mc10,
           sv_dir = "pcas/pca_referenceLayers", layers_format = ".asc",
           sv_proj_dir = "pcas/pca_proj")

# Read the pca object (output from ntbox function)
f1 <- readRDS("pcas/pca_referenceLayers/pca_object20_05_07_18_01.rds")

# Summary 
f2 <- summary(f1)

# The scree plot 
png(filename = "screeplot_merra.png", width = 1200*1.3, height = 1200*1.3, res = 300)
plot(f2$importance[3,1:5]*100, xlab = "Principal component", 
     ylab = "Percentage of variance explained", ylim=c(0,100),
     type = "b",frame.plot = T, cex = 1.5)
points(f2$importance[2,1:5]*100, pch = 17, cex = 1.5)
lines(f2$importance[2,1:5]*100, lty = 2, lwd = 1.5)
legend(x = 3.5, y = 60, legend = c("Cumulative", "Non-cumulative"),
       lty = c(1,2),pch = c(21,17),bty = "n",cex=0.85,pt.bg = 'white')
dev.off()

# PCs used were pc:1, 2, 3, 4
dir.create("Model_calibration/PCs_M")

nums <- 1:4

file.copy(from = paste0("pcas/pca_referenceLayers/PC0", nums,".asc"), 
          to = paste0("Model_calibration/PCs_M/PC0", nums,".asc"))

dir.create("pcas/pca_proj")

file.copy(from = paste0("pcas/pca_proj/PC0", nums,".asc"), 
          to = paste0("Model_calibration/PCs_world/PC0", nums,".asc"))
#-------------------------------------------------------------------------------


# The next parts help in performing Model calibration, Final model production, 
# and model projections.
# 1) We will prepare the data the way it needs to be prepared (SWD format). This format
# does not use raster layers for calibration, instead we will use the values of such layers 
# extracted to csv files.
# 2) Model calibration will be done using the previously prepared data.
# 3) Final models will be constructed and projected to the World using the
# parameters selected during model calibration.

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


# The next part helps to run all the processes described below 
# 1) Filter those replicates that could
#    predict invasion and native occurrence
#    points (training and testing)
#    the results will be store in 
#    01_medians_by_replicates (binary and median maps of replicates) and 
#    01_medians_by_replicates_results (information
#    on features classes and regulation multiplier
#    for the selected replicates)
# 2) Estimate the medians and sum of thresholded maps 
#    that could predict the by each dataset (1-5)
# 3) Filter variables that will be used to compute the mop (03_mops_by_dataset)
# 4) Estimate the medians of medians and the global sum of binary maps 

#--------------------------------Model consensus--------------------------------
# Data to test final models
## Invasive points
invs <- read.csv("coords_invasion_nortam.csv")

## All points of calibration area 
nats <- read.csv(list.files(dtrv, pattern = "all.csv$")[1])


# Consensus for distance thinned records raw variables
## Folder path where final models are
dtrv <- "Model_calibration/DT_RV"

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = dtrv)


# Consensus for distance thinned records raw variables
## Folder path where final models are
dtpc <- "Model_calibration/DT_PC"

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = dtpc)


# Consensus for distance thinned records raw variables
## Folder path where final models are
ctrv <- "Model_calibration/CT_RV"

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = ctrv)


# Consensus for distance thinned records raw variables
## Folder path where final models are
ctpc <- "Model_calibration/CT_PC"

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = ctpc)
#-------------------------------------------------------------------------------


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


#------------------------Preparing data for simulations-------------------------
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
