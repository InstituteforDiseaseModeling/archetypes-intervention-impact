## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 00_extraction_functions.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# May 2018
# 
# Collection of functions for covariate extraction and matrix manipulation. 
# 
# For a detailed project write-up see
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(data.table)
library(stringr)

get_mask <- function(continent, in_fname, out_fname){
  print("generating mask")
  mask_vals <- list('africa'=1,
                    'americas'=2,
                    'asia'=3)
  mask <- raster(in_fname)
  mask[mask!=mask_vals[[continent]]] <- NA
  clipped_mask <- trim(mask)
  writeRaster(clipped_mask, out_fname, overwrite=T)
  return(clipped_mask)
}

save_raster <- function(full_raster, out_fname, mask){
  # crop, save, and return values for a Raster* object
  vals <- extend(crop(full_raster, mask), mask)
  compareRaster(vals, mask)
  vals <- raster::mask(vals, mask, maskvalue=0)

  print("saving raster")
  writeRaster(vals, out_fname, bylayer=T, overwrite=T)
  
  return(vals)
}

# function to align resolutions between two rasters
align_res <- function(rast, template.rast){
  if (!identical(res(rast), res(template.rast))) {
    rast  <- raster::resample(rast, template.rast, method = 'ngb')
  }
  return(rast)
}



