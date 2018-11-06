## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# extract_and_svd.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# May 2018
# 
# For a given covariate and continent of interest, this script extracts covariate-specific
# data from global rasters, transforms it into a panel dataset, and saves it for further clustering analyses.  
# 
# For a detailed project write-up see
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------


library(gdistance)
library(data.table)
library(stringr)
library(stats)

rm(list=ls())

source("classify_functions.r")
root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")
mask_dir <- "Z:/mastergrids/Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif"

extra_crop_rasters <- list(asia="Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/5k/Synoptic/CHIRPS.Synoptic.01.mean.5km.Data.tif",
                           africa=file.path(base_dir, "vectors", "gambiae.tif")
                           )

overwrite <- F
cov_details <- fread("clustering_covariates.csv")

for (idx in 1:nrow(cov_details)){
  
  this_cov <- cov_details[idx]
  continents <- strsplit(this_cov$continents, "/")[[1]]
  
  for(continent in continents){
    
    print(paste("running extraction for", this_cov$cov, "in", continent))
    
    main_dir <- file.path(base_dir, continent, "rasters")
    dir.create(main_dir, showWarnings=F, recursive=T)
    all_vals_fname <- file.path(main_dir, "..", paste0(this_cov$cov, "_vals.csv"))
    
    ## find values
    if (file.exists(all_vals_fname) & overwrite==F){
      print("values already extracted")
    }else{
      
      print("clipping mask")
      if (continent %in% names(extra_crop_rasters)){
        mask <- get_mask(continent, out_dir=main_dir, mask_dir, extra_crop_rasters[[continent]])
      }else{
        mask <- get_mask(continent, out_dir=main_dir, mask_dir)
      }
      
      print("extracting from raster")
      all_vals <- lapply(strsplit(this_cov$values, "/")[[1]], extract_by_pattern, main_dir, this_cov, mask)
      all_vals <- rbindlist(all_vals)
      
      print("saving extracted values")
      write.csv(all_vals, file=all_vals_fname, row.names=F)
    }
    
  }

}

