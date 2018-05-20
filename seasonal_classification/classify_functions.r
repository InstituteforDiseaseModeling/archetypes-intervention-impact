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
  
  in_fname = paste0("TSI-Martens2-Pf.Synoptic.", str_pad(toString(month), 2, pad=0), ".Max.5km.Data.tif")
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