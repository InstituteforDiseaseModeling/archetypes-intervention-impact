library(gdistance)
library(data.table)
library(stringr)
library(stats)

rm(list=ls())

main_dir <- file.path(Sys.getenv("USERPROFILE"), "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/tsi_classification")
mask_dir <- "Z:/mastergrids/Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif"
tsi_dir <- "Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/"

all_tsi_dir <- file.path(main_dir, "tsi_vals.csv")

if (!file.exists(all_tsi_dir)){
  
  all_tsi <- lapply(1:12, function(month){
    
    print(paste("loading tsi for month", month))
    
    in_fname = paste0("TSI-Martens2-Pf.Synoptic.", str_pad(toString(month), 2, pad=0), ".Max.5km.Data.tif")
    out_fname = file.path(main_dir, "rasters", paste0("tsi_month_", month, ".tif"))
    
    if (file.exists(out_fname)){
      tsi <- raster(out_fname)
    }else{
      print("tsi not yet extracted! clipping global raster")
      
      if (!exists("afr_mask")){
        mask <- raster(mask_dir)
        afr <- mask == 1
        afr_mask <- trim(afr, values=F)
        rm(mask, afr)
      }
      
      tsi <- raster(file.path(tsi_dir, in_fname))
      tsi <- crop(tsi, afr_mask)
      writeRaster(tsi, out_fname)
    }
    tsi <- as.matrix(tsi)
    vals <- data.table(month=month, 
                       id = which(!is.na(tsi)),
                       tsi = tsi[!is.na(tsi)])

    return(vals)
  })
  
  all_tsi <- rbindlist(all_tsi)
  write.csv(all_tsi, file=all_tsi_dir, row.names=F)
  
}else{
  all_tsi <- fread(all_tsi_dir)
}

# format for svd
for_svd <- dcast(all_tsi, month ~ id)
svd_out <- svd(for_svd[, 2:ncol(for_svd)])
save(svd_out, file=file.path(main_dir, "svd_out.rdata"))

plot(svd_out$d^2/sum(svd_out$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "singular vectors", 
     ylab = "variance explained")

## keep only first n singular vectors
nvecs <- 2
sing_vecs <- svd_out$u[, 1:nvecs]

## multiply by original matrix to get rotations
rotation <- data.frame(t(t(sing_vecs)%*%as.matrix(for_svd[,2:ncol(for_svd)])))
rotation$id <- as.integer(rownames(rotation))
rotation <- data.table(rotation)
write.csv(rotation, file.path(main_dir, "svd_rotations.csv"), row.names=F)








