## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# extract_raster_values.r
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
library(rasterVis)

rm(list=ls())

source("classify_functions.r")

root_dir <- Sys.getenv("HOME")
map_root_dir <- "/Volumes/map_data/mastergrids"

if (!dir.exists(map_root_dir)){
  stop("Root map directory does not exist-- have you remembered to map the appropriate drives to your machine?")
}

# rewrite if there's already a saved covariate extraction?
overwrite_extraction <- F

base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/covariate_extraction")
mask_in_dir <- file.path(map_root_dir, "Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif")
transmission_limit_dir <- file.path(map_root_dir, "../GBD2017/Processing/Spatial_Data/Static_Limits", "Pf_Limits_EnvironmentalOnly_Endemic2017Only_5k.tif")

cov_details <- fread("clustering_covariates.csv")

# for now: don't mask based on transmission limits. Come back to this after replicating megatrends. 
# if (!file.exists(final_mask_dir) | overwrite_mask==T){
#   print("clipping full mask to transmission limits")
#   mask_layer <- raster(mask_in_dir)
#   mask_layer[mask_layer==0] <- -Inf
#   NAvalue(mask_layer) <- -Inf
#   writeRaster(mask_layer, temp_mask_dir, overwrite=T)
#   mask_layer <- raster(temp_mask_dir)
#   trans_limits <- raster(transmission_limit_dir)
#   mask_layer <- mask(mask_layer, trans_limits, maskvalue=0)
#   writeRaster(mask_layer, final_mask_dir, overwrite=T)
# }


for (idx in 1:nrow(cov_details)){
  
  this_cov <- cov_details[idx]
  continents <- strsplit(this_cov$continents, "/")[[1]]

  # append appropriate root dir to directories
  this_cov$dir <- file.path(map_root_dir, this_cov$dir)
  
  for(continent in continents){
    
    print(paste("running extraction for", this_cov$cov, "in", continent))
    
    main_dir <- file.path(base_dir, continent, "rasters")
    dir.create(main_dir, showWarnings=F, recursive=T)
    extraction_fname <- file.path(base_dir, continent, paste0(this_cov$cov, "_vals.csv"))
    
    ## find values
    if (file.exists(extraction_fname) & overwrite_extraction==F){
      print("values already extracted")
    }else{
      
      # check for mask
      clipped_mask_fname <- file.path(main_dir, "mask.tif")
      if (file.exists(clipped_mask_fname)){
        print("loading saved mask")
        mask <- raster(clipped_mask_fname)
      }else{
        print("clipping mask")
        mask <- get_mask(continent, in_fname=mask_in_dir, out_fname=clipped_mask_fname)
        plot(mask)
      }
      
      print("extracting from raster")
      all_vals <- lapply(strsplit(this_cov$values, "/")[[1]], extract_by_pattern, main_dir, this_cov, mask, overwrite_extraction)
      all_vals <- rbindlist(all_vals)
      
      print("saving extracted values")
      write.csv(all_vals, file=extraction_fname, row.names=F)
    }
    
  }

}

