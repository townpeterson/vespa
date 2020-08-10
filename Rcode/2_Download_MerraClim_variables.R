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

# Download bioclimatic variables

# needed package
library(raster)
library(rgdal)
library(rgeos)

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
# We need to set sea data as NA (the world shape file used here was: ESRI UIA_World Countries Boundaries, 
#2020. https://www.arcgis.com/home/item.html?id=252471276c9941729543be8789e06e12)
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

