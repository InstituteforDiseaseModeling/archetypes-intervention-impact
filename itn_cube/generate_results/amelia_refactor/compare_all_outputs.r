
library(data.table)
library(raster)


input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190517_refactor_database/"
func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
compare_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190507_sam_withseeds/"

source(file.path(func_dir, "check_file_similarity.r"))


## 01: check database
new_db <- fread(file.path(input_dir, "01_database.csv"))
old_db <- fread(file.path(compare_dir, "01_database.csv"))

old_db[, Tmean:=NULL]

check_sameness(old_db, new_db)

## 02: check covariates

new_covs <- fread(file.path(input_dir, "02_data_covariates.csv"))
load(file.path(compare_dir, "02_covariates.Rdata"))

old_covs <- data.table(all.covs)

check_sameness(old_covs, new_covs)


## 03: Access deviation



## 04: Use gap




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
