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

get_mask <- function(continent, main_dir, mask_dir){
  out_fname <- file.path(main_dir, "mask.tif")
  if (file.exists(out_fname)){
    clipped_mask <- raster(out_fname)
  }else{
    mask_vals <- list('africa'=1,
                      'americas'=2,
                      'asia'=3)
    mask <- raster(mask_dir)
    clipped_mask <- mask == mask_vals[[continent]]
    clipped_mask <- trim(clipped_mask, values=F)
    writeRaster(clipped_mask, out_fname)
  }
  return(clipped_mask)
}

extract_month <- function(month, main_dir, cov, mask_raster){
  
  print(paste("loading", names(cov), "for month", month))
  
  # find dataset to extract
  in_dir <- file.path(cov[[1]])
  files <- list.files(in_dir)
  monthval <-  str_pad(toString(month), 2, pad=0)
  pattern <- paste0(".*\\.", monthval, "\\..*.tif$")
  fname_idx <- which(str_detect(files, pattern))
  in_fname <- files[fname_idx]
  if (length(in_fname)>1){
    # try extracting mean value
    in_fname <- in_fname[in_fname %like% "mean"]
    
    if (length(in_fname)!=1){ # break if there are further issues
      stop(paste("multiple potential input rasters: ", in_fname))
    }
  }
  
  out_fname = file.path(main_dir, "rasters", paste0(names(cov), "_month_", month, ".tif"))
  
  if (file.exists(out_fname)){
    vals <- raster(out_fname)
  }else{
    print("clipping global raster")
    full <- raster(file.path(cov[[1]], in_fname))
    vals <- crop(full, mask_raster)
    mask_raster <- crop(mask_raster, vals) # ensure that extents of two rasters are the same
    vals <- mask(vals, mask_raster, maskvalue=FALSE)
    writeRaster(vals, out_fname)
  }
  vals <- as.matrix(vals)
  vals <- data.table(month=month, 
                     id = which(!is.na(vals)),
                     val = vals[!is.na(vals)])
  setnames(vals, "val", names(cov))
  
  return(vals)
  
}

rotate_matrix <- function(nvecs, main_dir, cov="tsi"){
  load(file.path(main_dir, paste0("svd_out_", cov, ".rdata")))
  sing_vecs <- svd_out$u[, 1:nvecs]
  
  ## multiply by original matrix to get rotations
  for_svd <- fread(file.path(main_dir, paste0(cov, "_vals.csv")))
  print("reshaping")
  for_svd <- dcast(for_svd, month ~ id)
  print("rotating")
  rotation <- data.frame(t(t(sing_vecs)%*%as.matrix(for_svd[,2:ncol(for_svd)])))
  rotation$id <- as.integer(rownames(rotation))
  rotation <- data.table(rotation)
  return(rotation)
}