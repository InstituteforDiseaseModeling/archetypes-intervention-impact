## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# joint_tsi_rainfall.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# June 2018
# 
# Extending the work described in the link below, this script takes the tsi and rainfall covariates 
# extracted in extract_ans_svd and runs svd on the joint dataset. 
# 
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------

library(gdistance)
library(data.table)
library(stringr)
library(stats)

rm(list=ls())

source("classify_functions.r")
base_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")
mask_dir <- "Z:/mastergrids/Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif"
cov_dirs <- list(tsi="Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/",
                 rainfall="Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/5k/Synoptic")
continents <- c("africa", "asia", "americas")

continents <- c("asia")

for (continent in continents){
  main_dir <- file.path(base_dir, continent)
  
  # load datasets, merge
  all_vals_fname <- file.path(main_dir, "joint_vals.csv")
  if (file.exists(all_vals_fname)){
    print("loading extracted values")
    all_vals <- fread(all_vals_fname)
  }else{
    tsi <- fread(file.path(main_dir, "tsi_vals.csv"))
    rainfall <- fread(file.path(main_dir, "rainfall_vals.csv"))
    
    # # create a joint sample raster with the correct extent
    # tsi_rast <- raster(file.path(main_dir, "rasters/tsi_month_1.tif"))
    # rainfall_rast <- raster(file.path(main_dir, "rasters/rainfall_month_1.tif"))
    # 
    all_vals <- merge(tsi, rainfall, by=c("month", "id"))
    all_vals <- melt(all_vals, id.vars=c("month", "id"), variable.name = "cov")
    write.csv(all_vals, file=all_vals_fname, row.names=F)
  }
  
  # svd
  svd_out_fname <- file.path(main_dir, "svd_out_joint.rdata")
  if (file.exists(svd_out_fname)){
    print("loading svd outputs")
    load(svd_out_fname)
  }else{
    print("running svd")
    for_svd <- dcast(all_vals, cov + month ~ id)
    svd_out <- svd(for_svd[, 3:ncol(for_svd)])
    save(svd_out, file=svd_out_fname)
  }
  
  png(file=file.path(main_dir, "figures", paste0("svd_out_joint_plot.png")))
  plot(svd_out$d^2/sum(svd_out$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "singular vectors",
       ylab = "variance explained")
  graphics.off()
  
}





