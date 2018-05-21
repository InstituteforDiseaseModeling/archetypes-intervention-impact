library(gdistance)
library(data.table)
library(stringr)
library(stats)

rm(list=ls())

base_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")
mask_dir <- "Z:/mastergrids/Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif"
cov_dirs <- list(tsi="Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/")
cov_dirs <- list(rainfall="Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/5k/Synoptic")
continents <- c("africa", "asia", "americas")
source("classify_functions.r")


for (idx in 1:length(cov_dirs)){
  
  this_cov <- cov_dirs[idx]
  
  for(continent in continents){
    
    print(paste("running extraction and svd  for", names(this_cov), "in", continent))
    
    main_dir <- file.path(base_dir, continent)
    all_vals_fname <- file.path(main_dir, paste0(names(this_cov), "_vals.csv"))
    
    ## find values
    if (file.exists(all_vals_fname)){
      print("loading extracted values")
      all_vals <- fread(all_vals_fname)
    }else{
      print("extracting from raster")
      mask <- get_mask(continent, file.path(main_dir, "rasters"), mask_dir)
      
      print("saving extracted values")
      all_vals <- lapply(1:12, extract_month, main_dir, this_cov, mask)
      all_vals <- rbindlist(all_vals)
      write.csv(all_vals, file=all_vals_fname, row.names=F)
    }
    
    # svd
    svd_out_fname <- file.path(main_dir, paste0("svd_out_", names(this_cov), ".rdata"))
    if (file.exists(svd_out_fname)){
      print("loading svd outputs")
      load(svd_out_fname)
    }else{
      print("running svd")
      for_svd <- dcast(all_vals, month ~ id)
      svd_out <- svd(for_svd[, 2:ncol(for_svd)])
      save(svd_out, file=svd_out_fname)
    }

    png(file=file.path(main_dir, "figures", paste0("svd_out_", names(this_cov), "_plot.png")))
      plot(svd_out$d^2/sum(svd_out$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "singular vectors",
           ylab = "variance explained")
    graphics.off()
    
  }

}

