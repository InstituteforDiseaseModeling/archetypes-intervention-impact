## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 01_extract_map_covariates.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# September 2019
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

source("extraction_functions.r")

root_dir <- Sys.getenv("HOME")
map_root_dir <- "/Volumes/map_data/mastergrids"

if (!dir.exists(map_root_dir)){
  stop("Root map directory does not exist-- have you remembered to map the appropriate drives to your machine?")
}

# rewrite if there's already a saved covariate extraction?
overwrite_extraction <- T

base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/covariates/no_transmission_limits")
mask_in_dir <- file.path(map_root_dir, "Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif")

cov_details <- fread("map_covariates.csv")

for (idx in 1:nrow(cov_details)){
  
  this_cov_details <- cov_details[idx]
  continents <- strsplit(this_cov_details$continents, "/")[[1]]

  # append appropriate root dir to directories
  this_cov_details$dir <- file.path(map_root_dir, this_cov_details$dir)
  
  for(continent in continents){
    
    print(paste("running extraction for", this_cov_details$cov, "in", continent))
    
    this_out_dir <- file.path(base_dir, continent, this_cov_details$cov)
    dir.create(this_out_dir, showWarnings=F, recursive=T)
    extraction_fname <- file.path(this_out_dir, paste0(this_cov_details$cov, "_vals.csv"))
    
    ## find values
    if (file.exists(extraction_fname) & overwrite_extraction==F){
      print("values already extracted")
    }else{
      
      # check for mask
      clipped_mask_fname <- file.path(base_dir, continent, "mask.tif")
      if (file.exists(clipped_mask_fname)){
        print("loading saved mask")
        mask_raster <- raster(clipped_mask_fname)
      }else{
        print("clipping mask")
        mask_raster <- get_mask(continent, in_fname=mask_in_dir, out_fname=clipped_mask_fname)
        plot(mask)
      }
      
      print("extracting from raster")
      
      all_vals <- lapply(unlist(strsplit(this_cov_details$variables, "/")), function(this_variable){
        
        print(paste("loading", this_cov_details$cov, "for", this_cov_details$variable_label, this_variable))
        
        # find raster to extract
        files <- list.files(this_cov_details$dir)
        pattern <- gsub("VAR", this_variable, this_cov_details$pattern)
        fname_idx <- which(str_detect(files, pattern))
        
        # make sure connection is strong enough to download file
        if (length(fname_idx)>1){ 
          stop(paste("multiple potential input rasters: ", this_variable))
        }else if (length(fname_idx)==0){
          print("raster initially not found, retrying")
          retries <- 10
          while(retries>0 & length(fname_idx)==0){
            print("retrying")
            Sys.sleep(1)
            files <- list.files(this_cov_details$dir)
            fname_idx <- which(str_detect(files, pattern))
            print(fname_idx)
            retries <- retries-1
          }
          if (length(fname_idx)==0){
            stop(paste("no input rasters found for pattern", pattern))
          }
        }
        
        # meat of the thing: cropping and extracting values from raster
        in_fname <- files[fname_idx]
        out_fname <- file.path(this_out_dir, paste0(this_cov_details$cov, "_", this_cov_details$variable_label, "_", this_variable, ".tif"))
        
        if (file.exists(out_fname) & overwrite_extraction==F){
          print("loading pre-extracted raster")
          val_raster <- raster(out_fname)
        }else{
          print("clipping global raster")
          val_raster <- save_raster(full_raster=raster(file.path(this_cov_details$dir, in_fname)), out_fname=out_fname, mask = mask_raster)
        }
      
        vals <- getValues(val_raster)
        vals <- data.table(cov = this_cov_details$cov,
                           variable_name = this_cov_details$variable_label,
                           variable_val = this_variable, 
                           id = which(!is.na(vals)),
                           value = vals[!is.na(vals)])
        
        return(vals)
      
      })
      
      all_vals <- rbindlist(all_vals)
      
      print("saving extracted values")
      write.csv(all_vals, file=extraction_fname, row.names=F)
    }
    
  }

}

