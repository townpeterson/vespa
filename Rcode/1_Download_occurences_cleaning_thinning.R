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

# Download GBIF data, data cleaning, thinning

#installing if needed
#install.packages("spocc")
#install.packages("rgbif")
#install.packages("maps")
#install.packages("rgdal")
#install.packages("raster")
#install.packages("sp")
#install.packages("dplyr")
#install.packages("purrr")
#
# url https://github.com/marlonecobos/ellipsenm#installing-the-package
#devtools::install_github("marlonecobos/ellipsenm")

# packages
library(spocc)
library(rgbif)
library(maps)
library(ellipsenm)
library(rgdal)
library(raster)
library(sp)
library(dplyr)
library(purrr)

# directory
setwd("Working_directory")

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

