###############################################################################################################
## 05_predict_itn.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt's original code to predict access and use outputs from the previously-run
## input models. 
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190508_replicate_inla func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data/ --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/05_predict_itn.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190508_replicate_inla/05_predictions --command 'Rscript ${CODE}'

rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190508_replicate_inla/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190508_replicate_inla/"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
  cov_dir <- "/Volumes/GoogleDrive/Team Drives/cubes/5km incomplete/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
  cov_dir <- Sys.getenv("cov_dir")
}

# directories
joint_dir <- file.path(joint_dir, "For_Prediction")

# TODO: load shapefiles and templates only when you're actually plotting (objects named World, Africa, Water, img)

# Load access and gap INLA outputs
# load(file.path(input_dir, "03_access_deviation.Rdata"))
# load(file.path(input_dir, "04_use_gap.Rdata"))

# load function script
source(file.path(func_dir, "05_predict_itn_functions.r"))
source(file.path(func_dir, "02_prep_covariate_functions.r"))

# why does *this* go to 2016 when our inla only goes to 2014/15?
prediction_years <- 2000:2016


## Load Covariates  TODO: Move all of this covariate extraction to 02_prep_covariates.r    ## ---------------------------------------------------------

cov_dt <- fread(file.path(func_dir, "covariate_key.csv"))
cov_dt[, used_sam:= as.logical(used_sam)]
cov_dt <- cov_dt[used_sam==T]

# find non-null raster indices-- we want to extract values for all of these every time
raster_indices <- which_non_null(file.path(joint_dir, '../african_cn5km_2013_no_disputes.tif'))


# Extract static covariates
static_fnames <- cov_dt[type=="static", list(fname=file.path(cov_dir, cov_name, fpath_append, fname))]

static_covs <- extract_values(static_fnames$fname, raster_indices)

# Extract annual covariates
covariate_years <- c(min=2001, max=2013) # this is the year range for which there are full monthly and annual time series

base_fnames <- lapply(cov_dt[type=="year"]$cov_name, get_annual_fnames, covariate_dt=cov_dt, input_dir=cov_dir)
base_fnames <- rbindlist(base_fnames)

# TODO: start year loop
this_year <- 2013

# make sure year is within the bounds of allowed covariate values
year_to_use <- min(this_year, covariate_years[["max"]]) 
year_to_use <- max(year_to_use, covariate_years[["min"]])

these_fnames <- copy(base_fnames)

these_fnames[, functional_year:=pmin(end_year, year_to_use)] # cap year by covariate availability
these_fnames[, functional_year:=pmax(functional_year, start_year)] 
these_fnames[, full_fname:=str_replace(base_fname, "YEAR", as.character(functional_year))]

annual_covs <- extract_values(these_fnames$full_fname, raster_indices, names=these_fnames$colname)

# Extract dynamic covariates

# TODO: month loop

this_month <- 3

these_fnames <- copy(cov_dt[type=="yearmon"])
these_fnames[, new_fname:=str_replace(fname, "YEAR", as.character(year_to_use))]
these_fnames[, new_fname:=str_replace(new_fname, "MONTH", str_pad(this_month, 2, pad="0"))]
these_fnames[, full_fname:=file.path(cov_dir, cov_name, fpath_append, new_fname)]

dynamic_covs <- extract_values(these_fnames$full_fname, raster_indices, names=these_fnames$cov_name)





