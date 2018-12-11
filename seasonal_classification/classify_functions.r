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

get_mask <- function(continent, in_fname, out_fname, extra_crop_dir=""){
  print("generating mask")
  mask_vals <- list('africa'=1,
                    'americas'=2,
                    'asia'=3)
  mask <- raster(in_fname)
  mask[mask!=mask_vals[[continent]]] <- NA
  clipped_mask <- trim(mask)
  if (extra_crop_dir!=""){
    extra_crop_raster <- raster(extra_crop_dir)
    extra_crop_raster <- crop(extra_crop_raster, clipped_mask)
    clipped_mask <- crop(clipped_mask, extra_crop_raster)
  }
  writeRaster(clipped_mask, out_fname, overwrite=T)
  return(clipped_mask)
}

extract_values <- function(raster_in_dir, raster_out_dir, mask, cov_name=""){
  full <- raster(raster_in_dir)
  vals <- crop(full, mask)
  if (!compareRaster(vals, mask, stopiffalse = F)){
    mask <- crop(mask, vals) # ensure that extents are the same. 
  }
  
  vals <- raster::mask(vals, mask, maskvalue=NA)
  # adjust null values in irs layer
  if (cov_name %in% c("irs_coverage", "act_coverage_older")){
    vals[vals==-9999] <- -Inf
  }
  # adjust values above 1 in act layer
  if (cov_name == "act_coverage"){
    vals[vals>1] <- 0
  }
  print("saving raster")
  writeRaster(vals, raster_out_dir, overwrite=T)
  
  return(vals)
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
    vals <- extract_values(raster_in_dir=file.path(cov$dir, in_fname), raster_out_dir=out_fname, mask = mask_raster, cov_name = cov$cov)
  }
  
  new_out_dir <- gsub("tif", "pdf", out_fname)
  print(paste("plotting raster to", new_out_dir))
  pdf(new_out_dir)
  print(plot(vals, main=paste(cov$cov, sweep_value)))
  graphics.off()
  
  plot(vals, main=paste(cov$cov, sweep_value))
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