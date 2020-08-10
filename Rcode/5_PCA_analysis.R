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

# PCA analysis 

## url https://github.com/luismurao/ntbox#complete-installation-guide
#remotes::install_github("luismurao/ntbox")
library(ntbox)

# stack of layers in M
varsm <- stack(list.files("calibration_area",pattern = ".asc$", full.names = T))

# stack of layers in the world
mc10 <- stack(list.files("merra_new", pattern = ".asc$", full.names = T))

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
#-------------------------------------------------------------------------------
