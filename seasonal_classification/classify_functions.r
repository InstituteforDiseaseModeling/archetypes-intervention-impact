## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# classify_functions.r
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

get_mask <- function(continent, out_dir, mask_dir, extra_crop_dir=NULL){
  out_fname <- file.path(out_dir, "mask.tif")
  if (file.exists(out_fname)){
    print("loading saved mask")
    clipped_mask <- raster(out_fname)
  }else{
    print("generating mask")
    mask_vals <- list('africa'=1,
                      'americas'=2,
                      'asia'=3)
    mask <- raster(mask_dir)
    clipped_mask <- mask == mask_vals[[continent]]
    clipped_mask <- trim(clipped_mask, values=F)
    if (!is.null(extra_crop_dir)){
      extra_crop_raster <- raster(extra_crop_dir)
      extra_crop_raster <- crop(extra_crop_raster, clipped_mask)
      clipped_mask <- crop(clipped_mask, extra_crop_raster)
    }
    writeRaster(clipped_mask, out_fname)
  }
  return(clipped_mask)
}

extract_values <- function(raster_in_dir, raster_out_dir, mask){
  full <- raster(raster_in_dir)
  vals <- crop(full, mask)
  # if (!compareRaster(vals, mask, stopiffalse = F)){
  #   mask <- crop(mask, vals) # ensure that extents are the same 
  # }
  vals <- raster::mask(vals, mask, maskvalue=FALSE)
  writeRaster(vals, raster_out_dir)
  return(vals)
}

extract_by_pattern <- function(sweep_value, out_dir, cov, mask_raster){
  
  print(paste("loading", cov$cov, "for", cov$variable, sweep_value))
  
  # find dataset to extract
  in_dir <- cov$dir
  files <- list.files(in_dir)
  pattern <- gsub("VAR", sweep_value, cov$pattern)
  fname_idx <- which(str_detect(files, pattern))
  in_fname <- files[fname_idx]
  if (length(in_fname)!=1){ # break if there are further issues
    stop(paste("multiple potential input rasters: ", in_fname))
  }

  out_fname <- file.path(out_dir, paste0(cov$cov, "_", cov$variable, "_", sweep_value, ".tif"))
  
  if (file.exists(out_fname)){
    vals <- raster(out_fname)
  }else{
    print("clipping global raster")
    vals <- extract_values(raster_in_dir=file.path(cov$dir, in_fname), raster_out_dir=out_fname, mask = mask_raster)
  }
  plot(vals)
  vals <- as.matrix(vals)
  vals <- data.table(cov = cov$cov,
                     variable_name = cov$variable,
                     variable_val = sweep_value, 
                     id = which(!is.na(vals)),
                     value = vals[!is.na(vals)])
  
  return(vals)
  
}


rotate_matrix <- function(nvecs, main_dir, cov="tsi"){
  load(file.path(main_dir, paste0(cov, "_svd.rdata")))
  sing_vecs <- svd_out$u[, 1:nvecs]
  
  ## multiply by original matrix to get rotations
  for_svd <- fread(file.path(main_dir, paste0(cov, "_vals.csv")))
  print("reshaping")
  for_svd <- dcast(for_svd, cov + variable_name + variable_val ~ id)
  print("rotating")
  rotation <- data.frame(t(t(sing_vecs)%*%as.matrix(for_svd[,4:ncol(for_svd)])))
  rotation$id <- as.integer(rownames(rotation))
  rotation <- data.table(rotation)
  return(rotation)
}