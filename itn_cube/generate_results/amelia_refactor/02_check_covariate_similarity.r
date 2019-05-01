
library(data.table)

rm(list=ls())

func_dir <-"/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/"
main_dir <- "/Users/bertozzivill/Desktop"

source(file.path(func_dir, "check_file_similarity.r"))

old <- load(file.path(main_dir, "preload.Rdata"))


new_static <- fread(file.path(main_dir, "static_covariates_27April2019.csv"))
new_annual <- fread(file.path(main_dir, "annual_covariates_27April2019.csv")) 
new_dynamic <- fread(file.path(main_dir, "dynamic_covariates_27April2019.csv"))

data <- data.table(data)

# compare static outputs
# cols_to_keep <- names(new_static)[!names(new_static) %in% names(data)]
# new_static <- new_static[, cols_to_keep, with=F]
# new_static[, fulldate:=NULL]

old_static <- data.table(static.covs)
old_static <- cbind(data, old_static)
old_static[, Survey:=as.character(Survey)]

check_sameness(old_static, new_static)

# compare annual outputs

old_annual <- data.table(yearonly.covs)
old_annual <- cbind(data, old_annual)
old_annual[, Survey:=as.character(Survey)]

new_names <- str_replace(names(new_annual), "IGBP_Landcover_", "landcover")
new_names <- str_replace(new_names, "AfriPop", "populatopn")
names(new_annual) <- new_names
new_annual <- new_annual[order(row_id)]

check_sameness(old_annual, new_annual)

# compare dynamic outputs
old_dynamic <- data.table(dynamic.covs)
old_dynamic <- cbind(data, old_dynamic)
old_dynamic[, Survey:=as.character(Survey)]

setnames(new_dynamic, c("EVI", "LST_day", "LST_night", "TCW", "TSI"), c("evy", "lst_day", "lst_night", "tcw", "tsi"))
new_dynamic <- new_dynamic[order(row_id)]

check_sameness(old_dynamic, new_dynamic)

