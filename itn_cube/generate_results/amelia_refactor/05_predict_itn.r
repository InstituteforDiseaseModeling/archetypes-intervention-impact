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


## Load Covariates  ## ---------------------------------------------------------

static_covs <- fread(file.path(input_dir, "02_static_covariates.csv"))
annual_covs <- fread(file.path(input_dir, "02_annual_covariates.csv"))
dynamic_covs <- fread(file.path(input_dir, "02_dynamic_covariates.csv"))






