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
continents <- c("africa", "asia", "americas")
cov_details <- fread("clustering_covariates.csv")

drop_na_cols = function(DT) {
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

for (continent in continents){
  main_dir <- file.path(base_dir, continent, "rasters")
  
  these_covs <- cov_details[continents %like% continent]
  full_label <- paste(these_covs$cov, collapse="_")
  
  # load datasets, merge
  all_vals_fname <- file.path(main_dir, paste0(full_label, ".csv"))
  if (file.exists(all_vals_fname)){
    print("loading extracted values")
    all_vals <- fread(all_vals_fname)
  }else{
    
    all_vals <- lapply(these_covs$cov, function(cov_name){
      vals <- fread(file.path(main_dir, paste0(cov_name, "_vals.csv")))
    })
    all_vals <- rbindlist(all_vals)
    write.csv(all_vals, file=all_vals_fname, row.names=F)
  }
  
  # svd
  svd_out_fname <- file.path(main_dir, "..", paste0(full_label, "_svd.rdata"))
  if (file.exists(svd_out_fname)){
    print("loading svd outputs")
    load(svd_out_fname)
  }else{
    print("running svd")
    for_svd <- dcast(all_vals, cov + variable_name + variable_val ~ id)
    drop_na_cols(for_svd) # drop columns with null values (due to extent differences in rasters)
    svd_out <- svd(for_svd[, 5:ncol(for_svd)])
    save(svd_out, file=svd_out_fname)
  }
  
  png(file=file.path(main_dir, paste0(full_label, "_svd.png")))
  plot(svd_out$d^2/sum(svd_out$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "singular vectors",
       ylab = "variance explained")
  graphics.off()
  
}



