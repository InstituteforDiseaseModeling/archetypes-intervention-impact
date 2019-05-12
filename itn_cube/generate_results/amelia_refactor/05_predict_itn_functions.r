###############################################################################################################
## 05_predict_itn_functions.r
## Amelia Bertozzi-Villa
## May 2019
## 
## Functions to accompany 05_predict_itn.r-- prepping survey data to go into the ITN cube model. 
## -which_non_null: way to find raster indices for which you want to extract other raster values
##############################################################################################################


which_non_null <- function(raster_fname){
  reference_raster <- raster(raster_fname)
  NAvalue(reference_raster) <- -9999
  reference_vals <- getValues(reference_raster)
  index_vals <- which(!is.na(reference_vals))
  return(index_vals)
}

