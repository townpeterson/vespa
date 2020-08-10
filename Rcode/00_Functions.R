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

# Necessary packages 
pkgs <- c("magrittr","raster","stringr","furrr","matrixStats")

# Check if they are installed
to_install <- pkgs[!pkgs %in% installed.packages()]

# Install missing packages
if(length(to_install)>0L) install.packages(pkgs)

# Loading packages 
library(magrittr)
library(raster)
library(stringr)
library(furrr)
library(matrixStats)
library(kuenm)


#-------------------------------------------------------------------------------
# Function to prepare swd format data
prep_swd <- function(vars, all, train, test, name_occ, back_out_dir) {
  # preparing swd format files
  ## adding variable values to records (changing coordinates to center of cells)
  ### all occurrences
  xyval <- extract(vars, all[, 2:3], cellnumbers = TRUE)
  xyras <- xyFromCell(vars, xyval[, 1])
  alld <- data.frame(sp=all[, 1], xyras, xyval[, -1])
  colnames(alld)[1:3] <- c("species", "longitude", "latitude")
  
  ### training set
  xyval <- extract(vars, train[, 2:3], cellnumbers = TRUE)
  xyras <- xyFromCell(vars, xyval[, 1])
  traind <- data.frame(sp=train[, 1], xyras, xyval[, -1])
  colnames(traind)[1:3] <- c("species", "longitude", "latitude")
  
  ### testing set
  xyval <- extract(vars, test[, 2:3], cellnumbers = TRUE)
  xyras <- xyFromCell(vars, xyval[, 1])
  testd <- data.frame(sp=test[, 1], xyras, xyval[, -1])
  colnames(testd)[1:3] <- c("species", "longitude", "latitude")
  
  ## saving records with variable values
  write.csv(alld, paste(name_occ, "_all.csv"), row.names = F)
  write.csv(traind, paste(name_occ, "_train.csv"), row.names = F)
  write.csv(testd, paste(name_occ, "_test.csv"), row.names = F)
  
  ## background files
  ### creating all combinations of variables
  var_names <- names(vars)
  var_sets <- all_var_comb(var.names = var_names, min.number = 2)
  
  ### adding all variable values to background
  backd <- data.frame(background = "background", rasterToPoints(vars))
  
  ### saving background files with distinct sets of variables
  bnames <- colnames(backd)[1:3] # three first columns
  sets <- names(var_sets) # distinct names of sets
  
  dir.create(back_out_dir) # new folder for all background
  
  #### writing all background files
  sv <- sapply(sets, function(x) {
    nms <- c(bnames, var_sets[[x]])
    write.csv(backd[, nms], file = paste0(back_out_dir, "/", x, ".csv"), 
              row.names = FALSE)
  })
  return(NULL)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Function to do calibration and projections together
maxent_calib_project <- function(all_occ, tra_occ, tes_occ, back_dir, 
                                 mx_path, batch_cal, candidate_models_dir,
                                 cal_res_dir, all_gvars, batch_fin, 
                                 final_models_dir) {
  # model calibration
  rm <- c(0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6) # regularization multipliers
  fc <- c("lq", "lp", "q", "qp", "lqp", "lqpt", "lqpth", "lqph") # feature classes

  calib <- kuenm_cal_swd(occ.joint = all_occ, occ.tra = tra_occ, 
                         occ.test = tes_occ, back.dir = back_dir, 
                         batch = batch_cal, out.dir.models = candidate_models_dir, 
                         reg.mult = rm, f.clas = fc, 
                         maxent.path = mx_path, kept = FALSE, out.dir.eval = cal_res_dir)
  
  
  # final model
  ## preparing variables for projections 
  ### selected sets
  sel <- read.csv(paste0(cal_res_dir, "/selected_models.csv"), stringsAsFactors = F)
  sets <- strsplit(sel[, 1], split = "_")
  
  ### new directory
  dir.create("G_variables")
  
  ### copying files in new directory
  gsets <- lapply(sets, function(z) {
    set <- paste(z[5:length(z)], collapse = "_")
    bset <- read.csv(paste0(back_dir, "/", set, ".csv")) 
    gvs <- colnames(bset)[-(1:3)]
    
    dir.create(paste0("G_variables/", set))
    inin <- paste0("G_variables/", set, "/projection")
    dir.create(inin)
    file.copy(from = paste0(all_gvars, "/", gvs, ".asc"),
              to = paste0(inin, "/", gvs, ".asc"))
    set
  })
  
  
  ## running final models 
  gv <- "G_variables" # sets of variables for projections

  fmod <- kuenm_mod_swd(occ.joint = all_occ, back.dir = back_dir, 
                        out.eval = cal_res_dir, batch = batch_fin, rep.n = 10, 
                        rep.type = "Bootstrap", jackknife = TRUE, 
                        out.format = "cloglog", project = TRUE, G.var.dir = gv,
                        ext.type = "all", maxent.path = mx_path, 
                        out.dir = final_models_dir)
}
#-------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# This function will compute the medians and sums of Final models
# for those replicates that can predict both
# dAll (training and testing data) and invasion points
#
# Parameters
# @param final_folder_path Path to the folder with final models
# @param invasion_df A data frame wit coordinates of invasion points 
#        two columns (longitude and latitude)
# @param all_data_df A data.frame with coordinates of native points 
#        three columns (Sp_name,longitude,latitude)
# @param save_dir Path to the directory where results will be saved

medians_by_set <- function(final_folder_path,
                           invasion_df,
                           all_data_df,
                           save_dir){
  
  
  invs <- invasion_df
  nats <- all_data_df
  
  # Paths to final models and its replicates
  finals_path <- list.files(final_folder_path,
                            pattern = "Final_models[0-9]$",
                            full.names = T)
  
  # List replicates of the final models
  
  paths_asc <- list.files(finals_path,"Vespa_mandarinia_[0-9]_projection.asc",
                          recursive = T,full.names = T)
  
  # Mask to perform calculations pixels with values 
  
  rmolde <- raster::raster(paths_asc[1])
  
  rmolde <- rmolde[]
  
  # Pixel ids with values
  
  id_nonas <- which(!is.na(rmolde))
  
  # Csv of the sample predictions
  
  paths_val <- list.files(finals_path,"Vespa_mandarinia_[0-9]_samplePredictions.csv",
                          recursive = T,full.names = T)
  # Joint paths to ascii files
  
  paths_levels <- paths_asc %>% 
    stringr::str_split(.,"/",simplify = T) %>% 
    data.frame(.,stringsAsFactors = F)
  level_names <- names(paths_levels)
  # Names of the models
  fmodels_cname <- level_names[grep(pattern = "Final*",
                                    x = paths_levels[1,])]
  
  models_cname <- level_names[grep(pattern = "^M_",
                                   x = paths_levels[1,])]
  
  dfLevels <- data.frame(paths_levels,paths_asc,
                         paths_val,
                         stringsAsFactors = F)
  # Split paths by thining method
  modelLevelsList <- split(dfLevels,dfLevels[,fmodels_cname])
  
  # Create a directory to save results
  
  if(!dir.exists(save_dir)) dir.create(save_dir)
  
  # Number final models 
  nfinals <- seq_along(modelLevelsList)
  
  # Run models in parallel
  plan(multisession)
  r1 <- nfinals %>% purrr::map_df(function(x){
    modelDF <- modelLevelsList[[x]]
    # Paths by extrapolation type
    modelDFL <- split(modelDF,modelDF[,models_cname])
    nextraSets <- seq_along(modelDFL)
    # Estimate model threshold at five percent 
    extraR <-   nextraSets %>% furrr::future_map_dfr(function(y){
      modelsetDF <- modelDFL[[y]]
      n_mod <- 1:nrow(modelsetDF)
      mods_test <- n_mod %>% purrr::map(function(z){
        rmod <- raster::raster(modelsetDF$paths_asc[z])
        xyVal <- read.csv(modelsetDF$paths_val[z])
        modpredVal <- sort(raster::extract(rmod,xyVal[,1:2]))
        nrecords <- nrow(xyVal) 
        fivepercent <- ceiling(nrecords*0.05)
        thFiveper <-  modpredVal[fivepercent]
        rbin <- rmod > thFiveper 
        
        # Check if models predict invasion points and native points at 5%

        binpreds <- raster::extract(rbin,nats[,2:3])
        invpreds <- raster::extract(rbin,invs[,1:2])
        #omr <- 1- sum(binpreds)/nrecords
        # Check which models could predict invasion points
        nopredicted <- which( binpreds ==0)
        if(all(invpreds==1) & (length(nopredicted) <=  fivepercent)){
          
          return(list(modelsetDF[z,],rbin))
        }
        
      })
      # List of replicates that could predict the invasion points
      dfmodsList <- lapply(seq_along(mods_test),
                           function(i) mods_test[[i]][[1]])
      binmodsList <- unlist(lapply(seq_along(mods_test),
                                   function(i) mods_test[[i]][[2]]))
      # Compute the medians for those replicates that could predict 
      # invasion points 
      if(length(binmodsList)>0L){
        # Stack if replicates that predicted invasion points
        bimdosStack <- raster::stack( binmodsList)
        # Sum  predicted cells
        valsbin <- rowSums(bimdosStack[id_nonas])
        # Data frame with replicate paths
        models_passDF <- do.call("rbind.data.frame",dfmodsList)
        models_bin_sum <- bimdosStack[[1]]
        models_bin_sum[id_nonas] <- valsbin
        models_cont_S <- raster::stack(models_passDF$paths_asc)
        # Compute medians
        medianVals <- matrixStats::rowMedians(models_cont_S[id_nonas])
        models_median <- models_cont_S[[1]]
        # Store medians in a raster
        models_median[id_nonas] <- medianVals
        # Paths of consensus maps 
        fname_base <- paste0(models_passDF[1,c(fmodels_cname,models_cname)]
                             ,collapse = "_")
        fname_base <- file.path(save_dir, fname_base)
        file_name_binsum <- paste0(fname_base,
                                   "_bin_sum_replicates.tif")
        file_name_median <- paste0(fname_base,
                                   "_median_replicates.tif")
        # Save results 
        raster::writeRaster(models_median,file_name_median,overwrite=T)
        raster::writeRaster(models_bin_sum,file_name_binsum,overwrite=T)
        
        return(models_passDF)
      }
      else
        return()
      
      
    },.progress = T)
    return(extraR)
  })
  plan(sequential)
  return(r1)
}
#-------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# The code computes the medians of medians for each data set
#
# Parameters
# @param medians_by_set Path to medians of models that could predict invasion 
#        points by dataset
# @param sv_mmedians Directory to save medians

medians_dataset <- function(medians_by_set,
                            sv_mmedians){
  tifs_medians <- list.files(medians_by_set,
                             pattern = "median_replicates.tif$",
                             full.names = T)
  
  #tifs_medians_path <- stringr::str_split(tifs_medians,"_",simplify = T)
  # Paths to binary sums of replicates
  
  tifs_bins <- list.files(medians_by_set,
                          pattern = "bin_sum_replicates.tif$",
                          full.names = T)
  
  nm <- stringr::str_split(tifs_bins,"_",simplify = T)
  # Model set ID
  id_mods <- grep("models",nm[1,])
  # datasets 
  datasets <- unlist(stringr::str_extract_all(nm [,id_mods],"[0-9]"))
  datasets <- unique(datasets)
  nsets <- seq_along(datasets)
  rmodel <-  raster::raster(tifs_bins[1])
  # Pixels with values
  id_nonas <- which(!is.na( rmodel[]))
  # Create a directory where results will be saved
  if(!dir.exists(sv_mmedians)) dir.create(sv_mmedians)
  # Run process in parallel
  plan(multisession)
  # Number of final models
  nsets <- seq_along(  nsets)
  resmeds <- nsets %>% furrr::future_map(function(x){
    # Paths to final models for each set
    set <- grep(tifs_medians,pattern = paste0("Final_models",datasets[x]),
                value = T)
    # Make a stack of models for each type of extrapolation
    try(setE <- raster::stack(grep(set,pattern = "_E_",value = T)),silent = T)
    try(setEC <- raster::stack(grep(set,pattern = "_EC_",value = T)),silent = T)
    try(setNE <- raster::stack(grep(set,pattern = "_NE_",value = T)),silent = T)
    # Path for extrolation results
    extra_path <-  file.path(sv_mmedians,paste0("median_",
                                                paste0("Final_models",x,
                                                       "_E",".tif")))
    
    # Path for extrapolation with clamp results
    
    exclam_path <-  file.path(sv_mmedians,paste0("median_",
                                                 paste0("Final_models",x,
                                                        "_EC",".tif")))
    # Path for no-extropolation results
    
    noext_path <-  file.path(sv_mmedians,paste0("median_",
                                                paste0("Final_models",x,
                                                       "_NE",".tif")))
    
    # Compute medians for Extrapolation
    try({
      msetE <- rmodel
      mediansEVal <- matrixStats::rowMedians(setE[id_nonas])
      msetE[id_nonas] <- mediansEVal
      raster::writeRaster(msetE, extra_path,overwrite=T)
    },silent = T)
    # Compute medians for extrapolation and clamping
    try({
      mediansECVal <- matrixStats::rowMedians(setEC[id_nonas])
      msetEC <- rmodel
      msetEC[id_nonas] <- mediansECVal
      raster::writeRaster(msetEC, exclam_path ,overwrite=T)
    },silent = T)
    # Compute medians for no extrapolation
    try({
      msetNE <- rmodel
      mediansNEVal <- matrixStats::rowMedians(setNE[id_nonas])
      msetNE[id_nonas] <- mediansNEVal
      raster::writeRaster(msetNE, noext_path ,overwrite=T)
    },silent = T)
    
    return()
  },.progress = T)
  
  # Binary models sums (consensus maps)
  resBinary <- nsets %>% furrr::future_map(function(x){
    
    # Paths to final models for each set
    set <- grep(tifs_bins,pattern = paste0("Final_models",datasets[x]),
                value = T)
    # Make a stack of models for each type of extrapolation
    try(setE <- raster::stack(grep(set,pattern = "_E_",value = T)),silent = T)
    try(setEC <- raster::stack(grep(set,pattern = "_EC_",value = T)),silent = T)
    try(setNE <- raster::stack(grep(set,pattern = "_NE_",value = T)),silent = T)
    
    # Path for extrapolation 
    
    extra_path <-  file.path(sv_mmedians,paste0("binarySum_",
                                                paste0("Final_models",x,
                                                       "_E",".tif")))
    
    # Path for extrapolation with clamp results
    
    exclam_path <-  file.path(sv_mmedians,paste0("binarySum_",
                                                 paste0("Final_models",x,
                                                        "_EC",".tif")))
    
    # Path for no extrapolation
    
    
    noext_path <-  file.path(sv_mmedians,paste0("binarySum_",
                                                paste0("Final_models",x,
                                                       "_NE",".tif")))
    # Binary sum for extrapolation
    try({
      msetE <- raster::calc(setE,sum)
      raster::writeRaster(msetE, extra_path,overwrite=T)
    })
    # Binary sum for extrapolation and clamp
    
    try({
      msetEC <- raster::calc(setEC,sum)
      raster::writeRaster(msetEC, exclam_path ,overwrite=T)
    })
    # Binary sum for no extrapolation
    
    try({
      msetNE <- raster::calc(setNE,sum)
      raster::writeRaster(msetNE, noext_path ,overwrite=T)
    })
    
    return()
  },.progress = T)
  
  plan(sequential)
  
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Function to check the lambdas files of
# those replicates that could predict our invasion points
# and then filter the variables that were used by maxent
#
# @param df_mods data.frame with the paths of those models that predicted
#                invasion points. The script will look for the lambdas file
#               of those models
# @param varnames Names of the layers that were used to model the potential
#                 distribution of the wasp.


mops_to_do <- function(df_mods, varnames=paste0("PC0",1:4)){
  #df_mods <- sapply(df_mods,as.character)
  #id_nm <- sapply(df_mods[1,], as.character)
  # Position ID in the df_mods objeect of the names of the models that predicted 
  # invasion points No extrapolation NE, Extrapolation E, Extrapolation_Clampling
  id_nm <- grep(pattern = ".asc", df_mods[1,])
  #df_mods[[id_nm]] <- as.character(df_mods[[id_nm]])
  # Lambdas files
  lams_nams <- gsub("_projection.asc",".lambdas",df_mods[[id_nm]])
  # Lambdas path 
  lams_path <- sapply(seq_along(lams_nams), function(x){
    paste0(paste0(df_mods[x,1:(id_nm-1)],collapse = "/"),"/",
           lams_nams[x])
  })
  # Position of the final models names in df_mods
  id_final <- grep(pattern = "Final_models",
                   df_mods[1,])
  
  # Levels to split the data (models) given final models folder,
  # and extrapolation method
  
  levs <- paste0(df_mods[[id_final]],"_",df_mods[["model_type"]])
  
  df_mods2 <- data.frame(lams_path,levs,stringsAsFactors = F)
  
  df_mods2L <-df_mods2 %>% split(.$levs)
  # This code looks for variables used in each model replicte
  
  varsdf <- seq_along(df_mods2L) %>% purrr::map_df(function(x){
    
    set_name <- names(df_mods2L[x])
    dfl <- df_mods2L[[x]]
    # Read the lambdas and look for the environmental layers used
    df_vars <- 1:nrow(dfl) %>% purrr::map_df(function(y){
      lambd <- readLines(as.character(dfl$lams_path[y]))
      garbID <- stringr::str_detect(lambd,"linearPredictorNormalizer")
      garbID <- which(garbID)
      lambd2 <- lambd[1:(garbID-1)]
      lambd3 <- stringr::str_split(lambd2,",",simplify = T)
      colnames(lambd3)<- c("Variables","lambda","var_min","var_max")
      df_lambs <- data.frame(lambd3,stringsAsFactors = F)
      df_lambs$lambda <- as.numeric(df_lambs$lambda)
      nonzeros <- which(df_lambs$lambda != 0)
      df_lambs <- df_lambs[nonzeros,]
      
    })
    
    vars_used <- unlist(sapply(varnames, function(x){
      varin <- grep(pattern = x,df_vars$Variables)
      if(length(varin)>0L) return(x)
    }))
    vars_used <- paste0(vars_used,collapse = ",")
    data.frame(set_name,vars_used,stringsAsFactors = F)
    
  })
  # Return data frame of lambdas and layers used 
  return(varsdf)
  #df_mods %>% 
  
}
#-------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function to compute the great median: We will compute the median
# across all data sets 
# @param medians_medians_path path to the medians of medians (medians
#                              of final models across the datasets 
#                              by extrapolation type)
# @param save_path path to a directory where results will be saved


great_median <- function(medians_medians_path,save_path){
  # List models 
  tifs <- list.files(medians_medians_path,pattern = ".tif$",full.names = T)
  # Binary models
  bins <- grep("bina",tifs,value = T)
  # Continuous models
  meds <- grep("median_Final",tifs,value = T)
  # Extrapolatio types
  extypes <- c("_E.tif$","_EC.tif$","_NE.tif$")
  
  # Create a directory for saving results
  
  if(!dir.exists(save_path)) dir.create(save_path)
  
  # Compute consensus models for binary outputs
  rbins <-  seq_along(extypes) %>% purrr::map(function(x){
    binsS <- raster::stack(grep(bins,pattern = extypes[x],value = T))
    binSum <- raster::calc(binsS,sum)
    fname <- file.path(save_path,paste0("binarySums_allSets",extypes[x],".tif"))
    fname <- gsub(".tif[$]","",fname)
    raster::writeRaster(binSum,fname,overwrite=T )
  })
  
  # Parallel process
  
  plan(multisession)
  # Compute great medians for each extrapolation type
  rmeds <- seq_along(extypes) %>% furrr::future_map(function(x){
    medsS <- raster::stack(grep(meds,pattern = extypes[x],value = T))
    mm <- raster::calc(medsS,median)
    fname <- file.path(save_path,paste0("median_allSets",extypes[x],".tif"))
    fname <-gsub(".tif[$]","",fname)
    raster::writeRaster(mm ,fname,overwrite=T )
  })
  plan(sequential)
  
}
#-------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function to process final models. 
# 1) It estimates models medians for replicates 
#    that could predict the invasion points
# 2) Compute the medians of median of those models that could predict 
#    invasion points for each of the five datasets
# 3) Check the lambdas files of those replicates that could predict our invasion
#    points and then filters the variables that were used by maxent. These variables
#    are then used to compute the MOPs
# 4) Computes the great median
# @param invas A data frame with coordinates of invasion points 
#        two columns (longitude and latitude)
# @param nats A data.frame with coordinates of native points 
#        three columns (Sp_name,longitude,latitude)
# @param finals_path Path to the folder with all final models folders

comp_medians <- function(invs,nats,finals_path){
  
  cat("
      
  #----------------------------------------------------------------------
  #                        *First part* 
  #----------------------------------------------------------------------
  # This code will compute the medians and sums of Final models
  # for those replicates that can predict both
  # dAll (training and testing data) and invasion points
  #----------------------------------------------------------------------
  
      
      ")
  
  
  
  # Folder where I will save the results 
  save_res <- "all_results"
  if(!dir.exists(save_res)) dir.create(save_res)
  save_dir <- file.path(save_res,"01_medians_by_replicates")
  if(!dir.exists(save_dir)) dir.create(save_dir)
  
  # The function that helps us to archive the first goal
  
  
  
  df_models <- medians_by_set(final_folder_path = finals_path,
                              invasion_df = invs,
                              all_data_df = nats,
                              save_dir = save_dir)
  
  
  # The following code will help to distinguish between
  # models with extrapolation, extrapolation, and clamping
  # and no extrapolation
  
  ndf <- names(df_models)
  
  nidex <- which(ndf %in% c("paths_asc","paths_val"))
  
  df_mods <- df_models[,-nidex]
  extrapolID <- grep("_E$",df_mods$X3)
  extraclamID <- grep("_EC$",df_mods$X3)
  noextraID <- grep("_NE$",df_mods$X3)
  df_mods$model_type <- NA
  df_mods$model_type[extrapolID] <- "Extrapolation"
  df_mods$model_type[extraclamID] <- "Extrapolation_Clampling"
  df_mods$model_type[noextraID] <- "No_Extrapolation"
  
  # Save the results of model replicates that could predict our 
  # occurence data points
  save_rep_res <- file.path(save_res,"01_medians_by_replicates_results")
  if(!dir.exists(save_rep_res)) dir.create(save_rep_res)
  selmod_path <- file.path(save_rep_res, 
                           "selected_model_replicates.csv")
  
  write.csv(df_mods,selmod_path,row.names = F)
  
  cat("
  #----------------------------------------------------------------------
  #                   *Second  part*
  # The code computes the medians of medians for each data set
  #----------------------------------------------------------------------
      
      ")
  
  
  medians_by_set <- save_dir
  sv_mmedians <-  file.path(save_res,"02_medians_by_dataset")
  if(!dir.exists(sv_mmedians)) dir.create(sv_mmedians )
  
  
  medians_of_meds <- medians_dataset(medians_by_set = medians_by_set,
                                     sv_mmedians = sv_mmedians)
  
  cat("
  #----------------------------------------------------------------------
  #                     *Third part* 
  #----------------------------------------------------------------------
  # In this part we will check the lambdas files of
  # those replicates that could predict our invasion points
  # and then filter the variables that were used by maxent
  #----------------------------------------------------------------------
  
      
      ")
  
  
  df_selected_mods <-read.csv(selmod_path,stringsAsFactors = F)
  
  
  modeling_vars <- names(nats)[-(1:3)]
  
  mopvariablesDF <- mops_to_do(df_mods = df_selected_mods,
                               varnames = modeling_vars)
  
  # Each element of the vector represent a MOP
  # that you should do 
  #mopvariables <- unique(mopvariablesDF$vars_used)
  mop_path <- file.path(save_res,"03_mops_by_dataset")
  if(!dir.exists(mop_path)) dir.create(mop_path)
  mop_file <- file.path(mop_path,"mops_table.csv")
  write.csv(mopvariablesDF,mop_file ,row.names = F)
  
  cat("
  #----------------------------------------------------------------------
  #                    *Fourth part* 
  #----------------------------------------------------------------------
  # The great median: We will compute the median
  # across all data sets (remember there are 5).
      
      ")
  
  
  
  medians_medians_path <-  sv_mmedians
  save_great_medians <- file.path(save_res,"04_medias_of_medians")
  
  rgreat_medians <- great_median(medians_medians_path = medians_medians_path,
                                 save_path = save_great_medians)
  
  cat("Check all_results folder!!!\n")
  cat("Folder *01_medians_by_replicates* are the medians of final\n" ,
      "models filtered by those that could predict native area and invasion points\n",
      "Check *01_medians_by_replicates_results* folder\n")
  
  cat("Folder *02_medians_by_replicates* are the medians of final\n" ,
      "models filtered by those that could predict native area and invasion points\n")
  
  cat("Folder *03_mop_by_dataset* has a .csv with mop variables to be\n",
      "used for each dataset" ,"\n")
  
  cat("Folder *04_medians_of_medians* has the medians of medians \n",
      "and the sums of binary models\n")
  
}



# Function to identify set for performing MOPs
id_mop_set <- function(back_dir, mop_table) {
  vsets <- unique(mop_table$vars_used)
  
  allbacks <- list.files(back_dir, pattern = ".csv$", full.names = T)
  abnames <- gsub(".csv$", "", list.files(back_dir, pattern = ".csv$"))
  
  vinsets <- lapply(1:length(abnames), function(x) {
    cn <- colnames(read.table(allbacks[x], header = T, sep = ",", nrows = 2))[, -(1:3)]
    paste0(cn, collapse = ",")
  })
  
  gsets <- sapply(vsets, function(x) {
    abnames[which(vinsets == x)]
  })
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
### Function used to select the suitability value that corresponds to the selected
# percentage threshold to be used in the binarization
# - median_mod - should be the output of an SDM
# - suits - are the suitability values assign by the model to each presence point
# - percent - should be a number between 0 and 100 that indicates the percentage of
#   presence points to be used for the binarization
### This function returns the suitability threshold to make the binarization of the
### continuous output map of the SDM


thfunc <- function(median_mod,suits,percent){
  # sort suitability values into ascending order and get rid of NAs
  suits <- na.omit(sort(suits))
  npts <- length(suits)
  # case 1:  all the suitability values are considered
  if(percent==0){
    threshold <- suits[1]
    #bin_mod <- median_mod>=threshold
    return(threshold)
  }
  # case 2: only a percentage of the suitability values are considered
  # define the smallest integer not less than the corresponding percentage of presence points
  pth <- ceiling(npts*percent)
  threshold <- suits[pth]
  
  #bin_mod <- median_mod>threshold
  return(threshold)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
### Function used to run simulations of the potential invasion paths of the species,
#   given an SDM and different values for the number of neighbor cells that the
#   species can visit in a time unit.
# - adL - A list with the adjacency matrices for different dispersal scenarios;
#         matrices are computed using the adj_mat function of the bam package.
#         The function requeires a raster map representing the SDM prediction
#         in a binary format and a numerical value indicating the neighbor cells
#         that the species is able to reach in a time unit.
# - bin_thmodels - binarized SDM outputs timed by the MOPs
# - base_name - character string that gives the name to the animation (gif) of  
#              the simulation
# - pasos - number of steps or time units to be simulated
### This function returns:
### (1) a list of binary maps each one representing the simulated, potential 
### dispersion of the species for a particular value of adL
### (2) an animation of the potential dynamics of invasion obteined from the
### previous list of maps


sim_disperal <- function(adL,bin_thmodels,base_name,pasos){
  # Run the process in parallel
  plan(multiprocess)
  # Increase the available memory to run process
  options(future.globals.maxSize = 3000 * 1024^2)
  # Run simulations for each dispersal scenario
  results_sims <- seq_along(adL) %>% furrr::future_map(function(x){
    # Number of theresholded maps trimed by the MOP (10)
    nlay <- raster::nlayers(bin_thmodels)
    # Adjacency matrix
    ad_matrix <- adL[[x]]
    # apply bam simulation process using occurrence points as the initial conditions
    sim_bin_thmodels <- lapply(1:nlay, function(y){
      sparse_mod <- bam::model2sparse(bin_thmodels[[y]])
      # Run the simulation 
      sim <- bam::sdm_sim(set_A = sparse_mod,
                          set_M =  ad_matrix,
                          initial_points = occs_sparse,
                          nsteps = pasos)
      # convert simulation results to raster
      sim_final <- bam::sim2Raster(sim,pasos)
      
    })
    sim_bin_thmodels <- raster::stack(sim_bin_thmodels)
    # Estimate concensus cells (cells that were
    # invaded in the dispersal scenarios)
    sim_bin_thmodels_sum <- raster::calc(sim_bin_thmodels,sum)
    
    return(sim_bin_thmodels_sum)
  },.progress = TRUE)
  
  # Create an animation with the maps obtained with the different distance values
  plan(sequential)
  names(results_sims) <- names(adL)
  if(!dir.exists(base_name)) dir.create(base_name)
  
  base_path <- normalizePath(base_name)
  base_path <- gsub("[\\]","/",base_path)
  # Animation path 
  anipath <- file.path( base_path,
                        paste0("animation_sims_",
                               base_name,".gif"))
  
  #Pathos to the raster files 
  rasters_path <- file.path( base_path,
                             paste0("Sim_",
                                    base_name,"_D_",
                                    names(adL),".tif"))
  # Save the animations 
  
  animation::saveGIF({
    for(i in seq_along(results_sims)){
      plot(results_sims[[i]], 
           main= paste("Dispersal distance ",names(results_sims[i])))
    }
  },movie.name = anipath,interval = 0.8,
  ani.width = 1200, ani.height = 1200, ani.res = 300)
  
  res_d <- lapply(seq_along(rasters_path), function(x){
    writeRaster(results_sims[[x]],  rasters_path[x],overwrite=T)
  })
  return(res_d)
}
#-------------------------------------------------------------------------------
