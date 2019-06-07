
###############################################################################################################
## compare_new_format.r
## Amelia Bertozzi-Villa
## June 2019
## 
## Test new and old results (but with the same output format) for equivalence
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-16 --logging gs://map_data_z/users/amelia/logs --input-recursive new_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor old_dir=gs://map_data_z/users/amelia/itn_cube/results/20190521_replicate_prediction --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/compare_new_format.r --output compare_out_path=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam/05_predictions/compare_newformat_tifs.pdf --command 'Rscript ${CODE}'

rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c( "raster", "data.table", "rasterVis", "stats", "RColorBrewer", "gridExtra"))

if(Sys.getenv("func_dir")=="") {
  new_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam/"
  old_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_rerun_sam/"
  compare_out_path <- file.path(new_dir, "05_predictions/compare_tifs.pdf")
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  new_dir <- Sys.getenv("new_dir")
  old_dir <- Sys.getenv("old_dir")
  compare_out_path <- Sys.getenv("compare_out_path") 
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}


source(file.path(func_dir, "check_file_similarity.r"))

## 01: check database
print("comparing old and new data files")
new_db <- fread(file.path(new_dir, "01_database.csv"))
old_db <- fread(file.path(old_dir, "01_database.csv"))

check_sameness(old_db, new_db)

## 02: check covariates
print("comparing old and new data covariate files")
new_covs <- fread(file.path(new_dir, "02_data_covariates.csv"))
old_covs <- fread(file.path(old_dir, "02_data_covariates.csv"))
check_sameness(old_covs, new_covs)


## 03: Access deviation
print("comparing old and new access deviation files")
load(file.path(new_dir, "03_access_deviation.Rdata"))

new_acc_fixed <- mod_pred_acc$summary.fixed
new_acc_fixed$cov<-rownames(new_acc_fixed)
new_acc_fixed <- data.table(new_acc_fixed)
new_acc_fixed <- new_acc_fixed[order(cov)]

new_hyperpar <- mod_pred_acc$summary.hyperpar
new_hyperpar$metric<-rownames(new_hyperpar)
new_hyperpar <- data.table(new_hyperpar)

load(file.path(old_dir, "03_access_deviation.Rdata"))
old_acc_fixed <- mod_pred_acc$summary.fixed
old_acc_fixed$cov<-rownames(old_acc_fixed)
old_acc_fixed <- data.table(old_acc_fixed)
old_acc_fixed <- old_acc_fixed[order(cov)]

old_hyperpar <- mod_pred_acc$summary.hyperpar
old_hyperpar$metric<-rownames(old_hyperpar)
old_hyperpar <- data.table(old_hyperpar)

print("fixed effects")
check_sameness(old_acc_fixed, new_acc_fixed, sameness_cutoff = 1e-04)
print("hyperparameters")
check_sameness(old_hyperpar, new_hyperpar, sameness_cutoff = 1e-04)

## 04: Use gap
print("comparing old and new use gap files")
load(file.path(new_dir, "04_use_gap.Rdata"))

new_use_fixed <- mod_pred_use$summary.fixed
new_use_fixed$cov<-rownames(new_use_fixed)
new_use_fixed <- data.table(new_use_fixed)
new_use_fixed <- new_use_fixed[order(cov)]

new_hyperpar <- mod_pred_use$summary.hyperpar
new_hyperpar$metric<-rownames(new_hyperpar)
new_hyperpar <- data.table(new_hyperpar)

load(file.path(new_dir, "04_use_gap.Rdata"))
old_use_fixed <- mod_pred_use$summary.fixed
old_use_fixed$cov<-rownames(old_use_fixed)
old_use_fixed <- data.table(old_use_fixed)
old_use_fixed <- old_use_fixed[order(cov)]

old_hyperpar <- mod_pred_use$summary.hyperpar
old_hyperpar$metric<-rownames(old_hyperpar)
old_hyperpar <- data.table(old_hyperpar)
print("fixed effects")
check_sameness(old_use_fixed, new_use_fixed, sameness_cutoff = 1e-04)
print("hyperparameters")
check_sameness(old_hyperpar, new_hyperpar, sameness_cutoff = 1e-04)


# 05: .tifs
print("comparing old and new raster files")
old_raster_dir <- file.path(old_dir, "05_predictions")
new_raster_dir <- file.path(new_dir, "05_predictions")

compare_tifs <- function(old_tif, new_tif, name="", cutoff=0.001){
  
  print(name)
  
  # new_tif[old_tif==0] <- 0
  # old_tif[new_tif==0] <- 0
  diff <- old_tif-new_tif
  old_diff <- old_tif[abs(diff)>cutoff]
  
  if (length(old_diff)>0){
    print(paste(length(old_diff), "values above cutoff"))
    if (max(old_diff)==0){
      print("Discrepancy comes from zeros in old tiff, no issue")
    }else{
      print("Problem: see discrepancies below")
      print(summary(old_diff))
    }
  }else{
    print("No discrepancy")
  }
  
  stacked <- stack(old_tif, new_tif)
  names(stacked) <- c(paste("Old", name), paste("New", name))
  
  stackplot <- levelplot(stacked,
                         par.settings=rasterTheme(region=brewer.pal(8, "RdYlGn")),
                         xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F)
  
  
  return(stackplot)
}


pdf(compare_out_path, width=11, height=7)
for (this_year in 2000:2016){
  print(this_year)
  old_mean_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".MEAN.tif")))
  new_mean_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".MEAN.tif")))
  mean_stack <- compare_tifs(old_mean_tif, new_mean_tif, name="National Access") 
  
  old_dev_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".DEV.tif")))
  new_dev_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".DEV.tif")))
  dev_stack <- compare_tifs(old_dev_tif, new_dev_tif, name="Access Dev") 
  
  old_access_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".ACC.tif")))
  new_access_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".ACC.tif")))
  access_stack <- compare_tifs(old_access_tif, new_access_tif, name="Access")
  
  old_gap_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".GAP.tif")))
  new_gap_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".GAP.tif")))
  use_gap_stack <- compare_tifs(old_gap_tif, new_gap_tif, name="Use Gap")
  
  old_use_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".USE.tif")))
  new_use_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".USE.tif")))
  use_stack <- compare_tifs(old_use_tif, new_use_tif, name="Use") 
  
  if (file.exists(file.path(new_raster_dir, paste0("ITN_", this_year, ".RAKED_USE.tif")))){
    old_raked_use_tif <- raster(file.path(old_raster_dir, paste0("ITN_", this_year, ".RAKED_USE.tif")))
    new_raked_use_tif <- raster(file.path(new_raster_dir, paste0("ITN_", this_year, ".RAKED_USE.tif")))
    raked_use_stack <- compare_tifs(old_raked_use_tif, new_raked_use_tif, name="Raked Use")
  }else{
    raked_use_stack <- use_stack
  }
  
  
  full_stack <- list(mean_stack, dev_stack, access_stack, use_gap_stack, use_stack, raked_use_stack)
  full_plot <- grid.arrange(grobs=full_stack, ncol=3, top=paste(this_year))
  
}
graphics.off()





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
