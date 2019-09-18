## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 00_bound_by_transmission.r
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

source("00_extraction_functions.r")

root_dir <- Sys.getenv("HOME")
map_root_dir <- "/Volumes/map_data/mastergrids"

if (!dir.exists(map_root_dir)){
  stop("Root map directory does not exist-- have you remembered to map the appropriate drives to your machine?")
}

# rewrite if there's already a saved covariate extraction?
overwrite_extraction <- F

base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/00_covariate_extraction")
full_cov_dir <- file.path(base_dir, "no_transmission_limits")
bounded_cov_dir <- file.path(base_dir, "with_transmission_limits")

transmission_limit_dir <- file.path(map_root_dir, "../GBD2017/Processing/Spatial_Data/Static_Limits", "Pf_Limits_EnvironmentalOnly_Endemic2017Only_5k.tif")
transmission_limits <- raster(transmission_limit_dir)

cov_details <- fread("clustering_covariates.csv")

for (idx in 1:nrow(cov_details)){
  
  this_cov <- cov_details[idx]
  continents <- strsplit(this_cov$continents, "/")[[1]]
  
  for(continent in continents){
    
    print(paste("running extraction for", this_cov$cov, "in", continent))
    
    main_dir <- file.path(bounded_cov_dir, continent, "rasters")
    dir.create(main_dir, showWarnings=F, recursive=T)
    extraction_fname <- file.path(bounded_cov_dir, continent, paste0(this_cov$cov, "_vals.csv"))
    this_cov_dir <- file.path(full_cov_dir, continent, "rasters")
    
    ## find values
    if (file.exists(extraction_fname) & overwrite_extraction==F){
      print("values already extracted")
    }else{
      
      # check for mask
      clipped_mask_fname <- file.path(main_dir, "bounded_mask.tif")
      if (file.exists(clipped_mask_fname)){
        print("loading saved mask")
        mask <- raster(clipped_mask_fname)
      }else{
        print("finding mask")
        orig_mask <- raster(file.path(this_cov_dir, "mask.tif"))
        clipped_transmission <- extend(crop(transmission_limits, orig_mask), orig_mask)
        bounded_mask <- raster::mask(orig_mask, clipped_transmission, maskvalue=0, updatevalue=0)
        writeRaster(bounded_mask, clipped_mask_fname, overwrite=T)
        plot(bounded_mask)
      }
      
      print("extracting from raster")
      all_cov_fnames <- list.files(this_cov_dir)
      these_cov_fnames <- all_cov_fnames[all_cov_fnames %like% this_cov$cov & all_cov_fnames %like% "tif"]
      
      all_vals <- lapply(these_cov_fnames, function(this_fname){
        vals <- extract_values(raster_in_dir=file.path(this_cov_dir, this_fname), out_fname=file.path(main_dir, this_fname), mask = bounded_mask, cov_name = this_cov, plot=F)
        vals <- as.matrix(vals)
        
        # use the known format of the file to find the month/year/etc associated with this raster
        sweep_value <- gsub(".tif", "", unlist(strsplit(this_fname, "_")))
        sweep_value <- sweep_value[length(sweep_value)]
        vals <- data.table(cov = this_cov$cov,
                           variable_name = this_cov$variable,
                           variable_val = sweep_value, 
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

