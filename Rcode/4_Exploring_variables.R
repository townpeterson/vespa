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

# Exploring variables in environment

# packages
library(ellipsenm)
library(raster)

# working directory
setwd("Working_directory")

# reading layers
variables <- stack(list.files("merra_calibration", pattern = ".asc$", full.names = T))
variables <- variables[[c(1, 10:15, 2:9)]]

# reading records
records <- read.csv("V_mandarinia_50km.csv")


#-------------------Exploring variables in environmental space------------------
# exploring the data in environmental space
## only temperature variables
explore_espace(data = records, species = "scientificName", longitude = "longitude",
               latitude = "latitude", raster_layers = variables[[1:9]], save = T,
               name = "Temperature_variables.pdf")

## only precipitation variables
explore_espace(data = records, species = "scientificName", longitude = "longitude",
               latitude = "latitude", raster_layers = variables[[10:15]], save = T,
               name = "Precipitation_variables.pdf")

## highly correlated set
hcor <- c(1, 6, 8, 9, 10, 11, 12, 14, 15)
explore_espace(data = records, species = "scientificName", longitude = "longitude",
               latitude = "latitude", raster_layers = variables[[hcor]], 
               save = T, name = "High-cor_variables.pdf")

## non-correlated sets
ncor_sets <- list(c(3, 5, 7, 11, 12), c(3, 4, 5, 11, 12), c(2, 3, 10, 13), 
                  c(3, 4, 10, 13), c(3, 4, 12, 13), c(3, 7, 12, 13), 
                  c(3, 4, 13, 15), c(3, 7, 13, 15))

explore <- lapply(3:length(ncor_sets), function(x) {
  explore_espace(data = records, species = "scientificName", longitude = "longitude",
                 latitude = "latitude", raster_layers = variables[[ncor_sets[[x]]]], 
                 save = T, name = paste0("Non-cor_variables", x, ".pdf"))
})


# exploring variable correlation in one plot for all
jpeg("corrplot_merra.jpg", width = 120, height = 120, units = "mm", res = 600)
par(cex = 0.8)
vcor <- variable_correlation(variables, save = T, name = "correlation_merra", 
                             corrplot = T, magnify_to = 3)
dev.off()
