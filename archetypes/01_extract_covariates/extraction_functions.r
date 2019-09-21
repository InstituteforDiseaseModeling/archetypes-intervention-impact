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

extract_values <- function(full_raster, out_fname, mask){
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


extract_by_pattern <- function(sweep_value, out_dir, cov, mask_raster, overwrite=F){
  
  print(paste("loading", cov$cov, "for", cov$variable, sweep_value))
  
  # find dataset to extract
  in_dir <- cov$dir
  files <- list.files(in_dir)
  pattern <- gsub("VAR", sweep_value, cov$pattern)
  fname_idx <- which(str_detect(files, pattern))

  if (length(fname_idx)>1){ 
    stop(paste("multiple potential input rasters: ", in_fname))
  }else if (length(fname_idx)==0){
    print("raster initially not found, retrying")
    retries <- 10
    while(retries>0 & length(fname_idx)==0){
      print("retrying")
      Sys.sleep(1)
      files <- list.files(in_dir)
      fname_idx <- which(str_detect(files, pattern))
      print(fname_idx)
      retries <- retries-1
    }
    if (length(fname_idx)==0){
      stop(paste("no input rasters found for pattern", pattern))
    }
  }

  in_fname <- files[fname_idx]
  out_fname <- file.path(out_dir, paste0(cov$cov, "_", cov$variable, "_", sweep_value, ".tif"))
  
  if (file.exists(out_fname) & overwrite==F){
    print("loading pre-extracted raster")
    vals <- raster(out_fname)
  }else{
    print("clipping global raster")
    vals <- extract_values(full_raster=raster(file.path(cov$dir, in_fname)), out_fname=out_fname, mask = mask_raster)
  }
  
  plot_out_fname <- gsub("tif", "pdf", out_fname)
  print(paste("plotting raster to", plot_out_fname))
  pdf(plot_out_fname)
  print(plot(vals, main=paste(cov$cov, sweep_value)))
  graphics.off()
  
  vals <- as.matrix(vals)
  vals <- data.table(cov = cov$cov,
                     variable_name = cov$variable,
                     variable_val = sweep_value, 
                     id = which(!is.na(vals)),
                     value = vals[!is.na(vals)])
  
  return(vals)
  
}
