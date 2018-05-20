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
  if (length(fname_idx)>1){
    stop(paste("multiple potential input rasters: ", in_fname))
  }
  
  out_fname = file.path(main_dir, "rasters", paste0(names(cov), "_month_", month, ".tif"))
  
  if (file.exists(out_fname)){
    vals <- raster(out_fname)
  }else{
    print("clipping global raster")
    full <- raster(file.path(cov[[1]], in_fname))
    vals <- crop(full, mask_raster)
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