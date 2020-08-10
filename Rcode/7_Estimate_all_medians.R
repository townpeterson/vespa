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

# Script to perform model consensus results

# Load the R script with functions to compute the medians
source("R_code_V_mandarinia/00_Functions.R")

# Set as working directory the folder where your final models are
setwd("working directory")


# Run all the process 
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


#-----------------------Preparing model consensus results-----------------------
# Data to test final models
## Invasive points
invs <- read.csv("coords_invasion_nortam.csv")

# Consensus for distance thinned records raw variables
## Folder path where final models are
dtrv <- "Model_calibration/DT_RV"

## All points of calibration area 
nats <- read.csv(list.files(dtrv, pattern = "all.csv$")[1])

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = dtrv)


# Consensus for distance thinned records pcs
## Folder path where final models are
dtpc <- "Model_calibration/DT_PC"

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = dtpc)


# Consensus for country thinned records raw variables
## Folder path where final models are
ctrv <- "Model_calibration/CT_RV"

## All points of calibration area 
nats <- read.csv(list.files(ctrv, pattern = "all.csv$")[1])

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = ctrv)


# Consensus for country thinned records pcs
## Folder path where final models are
ctpc <- "Model_calibration/CT_PC"

# Median calculation
resMeds <- comp_medians(invs = invs, nats = nats, finals_path = ctpc)
#-------------------------------------------------------------------------------
