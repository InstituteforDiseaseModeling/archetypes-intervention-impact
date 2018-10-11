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
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification/")
continents <- c("africa", "asia", "americas")
continents <- c("africa")
cov_details <- fread("clustering_covariates.csv")
cov_details <- cov_details[variable=="month"]
overwrite <- T

for (continent in continents){
  print(paste("running svd on", continent))
  main_dir <- file.path(base_dir, continent)
  
  these_covs <- cov_details[continents %like% continent]
  full_label <- paste(these_covs$cov, collapse="_")
  
  # load datasets, merge
  all_vals_fname <- file.path(main_dir, paste0(full_label, "_vals.csv"))
  if (file.exists(all_vals_fname) & overwrite==F){
    print("loading extracted values")
    all_vals <- fread(all_vals_fname)
  }else{
    print("appending datasets")
    all_vals <- lapply(these_covs$cov, function(cov_name){
      vals <- fread(file.path(main_dir, paste0(cov_name, "_vals.csv")))
      # # temp
      # setnames(vals, c(cov_name, "month"), c("value", "variable_val"))
      # vals[, cov:=cov_name]
      # vals[, variable_name:="month"]
      return(vals)
    })
    
    # keep only those pixels with values for all covariates
    non_null_ids <- lapply(all_vals, function(df){
      return(unique(df$id))
    })
    shared_ids <- Reduce(intersect, non_null_ids)
    all_vals <- rbindlist(all_vals)
    all_vals <- all_vals[id %in% shared_ids]
    write.csv(all_vals, file=all_vals_fname, row.names=F)
  }
  
  # svd
  svd_out_fname <- file.path(main_dir, paste0(full_label, "_svd.rdata"))
  if (file.exists(svd_out_fname) & overwrite==F){
    print("loading svd outputs")
    load(svd_out_fname)
  }else{
    print("reshaping and filling nulls")
    for_svd <- dcast(all_vals, cov + variable_name + variable_val ~ id)
    print("running svd")
    svd_out <- svd(for_svd[, 4:ncol(for_svd)])
    save(svd_out, file=svd_out_fname)
  }
  
  print("plotting")
  png(file=file.path(main_dir, paste0(full_label, "_svd.png")))
  plot(svd_out$d^2/sum(svd_out$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "singular vectors",
       ylab = "variance explained")
  graphics.off()
  
}



