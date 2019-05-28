
library(data.table)
library(raster)
library(rasterVis)
library(stats)

rm(list=ls())

new_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190521_replicate_prediction/"
old_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190507_sam_withseeds/"
func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"

source(file.path(func_dir, "check_file_similarity.r"))

# load the .rdata objects that will disrupt the environment
load(file.path(old_dir, "03_access_deviation.Rdata"))
old_accdev <- copy(mod.pred)
load(file.path(old_dir, "04_use_gap.Rdata"))
old_use_gap <- copy(mod.pred)

## 01: check database
new_db <- fread(file.path(new_dir, "01_database.csv"))
old_db <- fread(file.path(old_dir, "01_database.csv"))

old_db[, Tmean:=NULL]

check_sameness(old_db, new_db)

## 02: check covariates

new_covs <- fread(file.path(new_dir, "02_data_covariates.csv"))
load(file.path(old_dir, "02_covariates.Rdata"))

old_covs <- data.table(all.covs)

check_sameness(old_covs, new_covs)


## 03: Access deviation
load(file.path(new_dir, "03_access_deviation.Rdata"))

new_acc_fixed <- mod_pred_acc$summary.fixed
new_acc_fixed$cov<-rownames(new_acc_fixed)
new_acc_fixed <- data.table(new_acc_fixed)
new_acc_fixed <- new_acc_fixed[order(cov)]

new_hyperpar <- mod_pred_acc$summary.hyperpar
new_hyperpar$metric<-rownames(new_hyperpar)
new_hyperpar <- data.table(new_hyperpar)

old_acc_fixed <- old_accdev$summary.fixed
old_acc_fixed$cov<-rownames(old_acc_fixed)
old_acc_fixed <- data.table(old_acc_fixed)
old_acc_fixed <- old_acc_fixed[order(cov)]

old_hyperpar <- old_accdev$summary.hyperpar
old_hyperpar$metric<-rownames(old_hyperpar)
old_hyperpar <- data.table(old_hyperpar)

check_sameness(old_acc_fixed, new_acc_fixed, sameness_cutoff = 1e-04)
check_sameness(old_hyperpar, new_hyperpar, sameness_cutoff = 1e-04)

## 04: Use gap
load(file.path(new_dir, "04_use_gap.Rdata"))

new_use_fixed <- mod_pred_acc$summary.fixed
new_use_fixed$cov<-rownames(new_use_fixed)
new_use_fixed <- data.table(new_use_fixed)
new_use_fixed <- new_use_fixed[order(cov)]

new_hyperpar <- mod_pred_acc$summary.hyperpar
new_hyperpar$metric<-rownames(new_hyperpar)
new_hyperpar <- data.table(new_hyperpar)

old_use_fixed <- old_accdev$summary.fixed
old_use_fixed$cov<-rownames(old_use_fixed)
old_use_fixed <- data.table(old_use_fixed)
old_use_fixed <- old_use_fixed[order(cov)]

old_hyperpar <- old_accdev$summary.hyperpar
old_hyperpar$metric<-rownames(old_hyperpar)
old_hyperpar <- data.table(old_hyperpar)

check_sameness(old_use_fixed, new_use_fixed, sameness_cutoff = 1e-04)
check_sameness(old_hyperpar, new_hyperpar, sameness_cutoff = 1e-04)


# 05: .tifs

this_year <- 2003

old_raster_dir <- file.path(old_dir, "05_predictions")
new_raster_dir <- file.path(new_dir, "05_predictions")

compare_tifs <- function(old_tif, new_tif, cutoff=0.001){
  
  new_tif[old_tif==0] <- 0
  old_tif[new_tif==0] <- 0
  diff <- old_tif-new_tif
  old_diff <- old_tif[abs(diff)>cutoff]
  
  if (length(old_diff)>0){
    print(paste(length(old_diff), "values above cutoff"))
    if (max(old_diff==0)){
      print("Discrepancy comes from zeros in old tiff, no issue")
    }else{
      print("Problem: see discrepancies below")
      print(summary(old_diff))
    }
  }else{
    print("No discrepancy")
  }
  
  print(levelplot(stack(old_tif, new_tif)))
  return(diff)
}


old_access_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".ACC.tif")))
new_access_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".ACC.tif")))
compare_tifs(old_access_tif, new_access_tif)


old_dev_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".DEV.tif")))
old_dev_tif <- calc(old_dev_tif, plogis)
new_dev_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".DEV.tif")))
dev_diff <- compare_tifs(old_dev_tif, new_dev_tif) # emplogit is weird so I can live with this

old_mean_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".MEAN.tif")))
new_mean_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".MEAN.tif")))
mean_diff <- compare_tifs(old_mean_tif, new_mean_tif) # this looks real, possibly related to how islands are treated

old_use_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".USE.tif")))
new_use_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, "USE")))
compare_tifs(old_use_tif, new_use_tif) # added instead of subtracting; fixing now

old_gap_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".GAP.tif")))
old_gap_tif <- calc(old_gap_tif, plogis)
new_gap_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".GAP.tif")))
levelplot(stack(old_gap_tif, new_gap_tif))
compare_tifs(old_gap_tif, new_gap_tif)



## Archival from covariate discrepancy:

# test <- data[evy!=compare_data$evy]
# compare_test <- compare_data[evy!=data$evy]
# 
# incorrect_points <- data.table(cellnumber=c(901138, 1027886, 1603132, 2141192),
#                                flooryear=c(2008, 2014, 2014, 2004),
#                                month=c(12, 2, 10, 12),
#                                sub_cellnumber=c(896098, 1026206, 1591372, 2137832))
# 
# dyn_data <- merge(compare_data, incorrect_points, by=c("cellnumber", "flooryear", "month"), all.x=T)
# dyn_data[is.na(sub_cellnumber), sub_cellnumber:=cellnumber]
