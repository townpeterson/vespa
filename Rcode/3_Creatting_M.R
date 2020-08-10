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

# Creating M based on Buffers and masking variables with M
library(raster)
library(rgeos)
library(rgdal)
library(maps)

# directory
setwd("Working_directory")

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
# Masking layers to M
list.files("calibration_area/")
# wait for the final M
# reading M again
M <- readOGR("calibration_area", layer = "M")

# if needed load variables form Merra

# masking layers to M
varsm <- mask(crop(mc10, M), M)

# saving masked layers as ascii files
lapply(names(varsm), function(x){
  writeRaster(varsm[[x]], paste0("calibration_area/",x,".asc"), overwrite=T)
})
#-------------------------------------------------------------------------------

