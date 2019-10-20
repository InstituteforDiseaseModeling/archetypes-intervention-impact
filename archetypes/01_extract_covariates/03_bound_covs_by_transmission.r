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

source("extraction_functions.r")

root_dir <- Sys.getenv("HOME")
# map_root_dir <- "/Volumes/map_data/mastergrids"
# 
# if (!dir.exists(map_root_dir)){
#   stop("Root map directory does not exist-- have you remembered to map the appropriate drives to your machine?")
# }

# rewrite if there's already a saved covariate extraction?
overwrite_extraction <- T

base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/covariates")
unbounded_cov_dir <- file.path(base_dir, "no_transmission_limits")
bounded_cov_dir <- file.path(base_dir, "with_transmission_limits")

# transmission_limit_dir <- file.path(map_root_dir, "../GBD2017/Processing/Spatial_Data/Static_Limits", "Pf_Limits_EnvironmentalOnly_Endemic2017Only_5k.tif")
transmission_limit_dir <- file.path(bounded_cov_dir, "from_GBD2017_Pf_Limits_EnvironmentalOnly_Endemic2017Only_5k.tif")
transmission_limits <- raster(transmission_limit_dir)

continents <- list.dirs(unbounded_cov_dir, recursive = F, full.names = F)

for (continent in continents){
  
  unbounded_continent_dir <- file.path(unbounded_cov_dir, continent)
  covs_to_bound <- list.dirs(unbounded_continent_dir, recursive = F, full.names = F)
  
  bounded_continent_dir <- file.path(bounded_cov_dir, continent)
  dir.create(bounded_continent_dir, showWarnings=F, recursive=T)
  
  for (this_cov in covs_to_bound){
    
    print(paste("running extraction for", this_cov, "in", continent))
    
    this_unbounded_cov_dir <- file.path(unbounded_continent_dir, this_cov)
    
    this_bounded_cov_dir <- file.path(bounded_continent_dir, this_cov)
    dir.create(this_bounded_cov_dir, showWarnings=F, recursive=T)
    
    unbounded_extraction_vals <- fread(file.path(this_unbounded_cov_dir, paste0(this_cov, "_vals.csv")))
    bounded_extraction_fname <- file.path(this_bounded_cov_dir, paste0(this_cov, "_vals.csv"))
    
    ## find values
    if (file.exists(bounded_extraction_fname) & overwrite_extraction==F){
      print("values already extracted")
    }else{
      
      # check for mask
      clipped_mask_fname <- file.path(bounded_continent_dir, "bounded_mask.tif")
      if (file.exists(clipped_mask_fname)){
        print("loading saved mask")
        bounded_mask <- raster(clipped_mask_fname)
      }else{
        print("finding mask")
        orig_mask <- raster(file.path(unbounded_continent_dir, "mask.tif"))
        clipped_transmission <- extend(crop(transmission_limits, orig_mask), orig_mask)
        bounded_mask <- raster::mask(orig_mask, clipped_transmission, maskvalue=0, updatevalue=0)
        writeRaster(bounded_mask, clipped_mask_fname, overwrite=T)
        plot(bounded_mask)
      }
      
      print("extracting from raster")
      all_cov_fnames <- list.files(this_unbounded_cov_dir)
      these_cov_fnames <- all_cov_fnames[all_cov_fnames %like% this_cov & all_cov_fnames %like% "tif"]
      
      this_cov_stack <- stack(file.path(this_unbounded_cov_dir, these_cov_fnames))
      cropped_layers <- crop_raster(this_cov_stack, mask = bounded_mask, 
                                    out_fname = file.path(this_bounded_cov_dir, paste0(names(this_cov_stack), ".tif"))
                                    )
      
      # convert to data table, find pixels to keep from originally extracted dataset.
      vals <- as.matrix(cropped_layers)
      vals <- data.table(vals,
                         id=1:nrow(vals))
      vals <- vals[complete.cases(vals)]
      vals <- unbounded_extraction_vals[id %in% unique(vals$id)]
      
      print("saving extracted values")
      write.csv(vals, file=bounded_extraction_fname, row.names=F)
    }
    
  }
  
}


