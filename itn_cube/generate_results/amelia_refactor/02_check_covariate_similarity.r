
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
cols_to_keep <- names(new_static)[!names(new_static) %in% names(data)]
new_static <- new_static[, cols_to_keep, with=F]
new_static[, fulldate:=NULL]

old_static <- data.table(static.covs)

check_sameness(old_static, new_static)

# compare annual outputs
cols_to_keep <- names(new_annual)[!names(new_annual) %in% names(data)]
new_annual <- new_annual[, cols_to_keep, with=F]
new_annual[, fulldate:=NULL]
names(new_annual) <- c("populatopn", "landcover0", "landcover1", "landcover2", "landcover3", "landcover4",
                       "landcover5", "landcover6", "landcover7", "landcover8", "landcover9", "landcover10",
                       "landcover11", "landcover12", "landcover13", "landcover14", "landcover15", "landcover16")
new_annual[, landcover13:=NULL] # removed from sam's dataset for population collinearity

old_annual <- data.table(yearonly.covs)

check_sameness(old_annual, new_annual)


# compare dynamic outputs
cols_to_keep <- names(new_dynamic)[!names(new_dynamic) %in% names(data)]
new_dynamic <- new_dynamic[, cols_to_keep, with=F]
new_dynamic[, fulldate:=NULL]
names(new_dynamic) <- c("evy", "lst_day", "lst_night", "tcw", "tsi")

old_dynamic <- data.table(dynamic.covs)

check_sameness(old_dynamic, new_dynamic)

