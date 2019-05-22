###############################################################################################################
## 05_predict_itn.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt's original code to predict access and use outputs from the previously-run
## input models. 
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190521_replicate_prediction func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data/ --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/05_predict_itn.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190521_replicate_prediction/05_predictions --command 'Rscript ${CODE}'

rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190517_refactor_database/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190517_refactor_database/"
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
set.seed(212)

# TODO: load shapefiles and templates only when you're actually plotting (objects named World, Africa, Water, img)

# Load access and gap INLA outputs
# load(file.path(input_dir, "03_access_deviation.Rdata"))
# load(file.path(input_dir, "04_use_gap.Rdata"))

# load function script
source(file.path(func_dir, "05_predict_itn_functions.r"))
source(file.path(func_dir, "01_create_database_functions.r")) # for calc_access
source(file.path(func_dir, "03_access_dev_functions.r")) # for 1-val emplogit

# why does *this* go to 2016 when our inla only goes to 2014/15?
prediction_years <- 2000:2016


## Load Covariates  ## ---------------------------------------------------------

this_year <- 2007

developing <- T

if (!developing){
  static_covs <- fread(file.path(input_dir, "02_static_covariates.csv"))
  annual_covs <- fread(file.path(input_dir, "02_annual_covariates.csv"))
  dynamic_covs <- fread(file.path(input_dir, paste0("02_dynamic_covariates/dynamic_", this_year, ".csv")))
  
  ## Get locations in x-y-z space of each pixel centroid for prediction ## ---------------------------------------------------------
  national_raster_dir <- file.path(joint_dir, "../african_cn5km_2013_no_disputes.tif")
  national_raster <- raster(national_raster_dir)
  
  prediction_indices <- which_non_null(national_raster_dir)
  prediction_cells <- data.table(row_id=prediction_indices, gaul=extract(national_raster, prediction_indices))
  
  prediction_cells <- cbind(prediction_cells, data.table(xyFromCell(national_raster, prediction_indices)))
  setnames(prediction_cells, c("x", "y"), c("longitude", "latitude"))
  prediction_xyz <- ll.to.xyz(prediction_cells)
}



## Get national access means from stock and flow  ## ---------------------------------------------------------

orig_stock_and_flow <- fread(file.path(input_dir, "01_stock_and_flow_prepped.csv"))
stock_and_flow <- orig_stock_and_flow[year>=this_year & year<(this_year+1)] # keep only the years we want to predict
iso_gaul_map <- fread(file.path(joint_dir, "National_Config_Data.csv"))
iso_gaul_map <- iso_gaul_map[ISO3 %in% stock_and_flow$iso3, list(ISO3, GAUL_Code, Name=IHME_Location_ASCII_Name)]

# load household size distributions for each survey
hh_sizes<-fread(file.path(joint_dir, 'Bucket_Model/HHsize.csv'))
survey_key=fread(file.path(joint_dir, 'Bucket_Model/KEY.csv'))

# format hh_sizes
hh_sizes[, V1:=NULL]
hh_sizes <- melt.data.table(hh_sizes, id.vars="HHSurvey", variable.name="hh_size", value.name="prop")
hh_sizes[, hh_size:=as.integer(hh_size)]

# BUG: The original script seems like it subsets on country-specific survey, 
# but I think it always uses the entire dataset to calculate home distributions because it filters on ISO instead of full country name
# (and also )

# update country names to fit with iso_gaul_map (for use when you actually DO calculate props by country name)
survey_key[, Status:=NULL]
setnames(survey_key, "Svy Name", "HHSurvey")
survey_key[Name=="Coted'Ivoire", Name:="Cote d'Ivoire"]
survey_key[Name=="Dem. Rep. of Congo", Name:="Democratic Republic of the Congo"]
survey_key[Name=="SaoTome & Principe", Name:="Sao Tome And Principe"]

# add two surveys not present in key
survey_key <- rbind(survey_key, data.table(HHSurvey=c("Kenya 2007", "Rwanda 2013"),
                                           Name=c("Kenya", "Rwanda")))

# merge with name map and survey distribution
survey_key <- merge(survey_key, iso_gaul_map, by="Name", all.x=T)
hh_sizes <- merge(hh_sizes, survey_key, by="HHSurvey", all.x=T)

# find size distribution across full dataset (make country-specific later);
# collapse such that the final bin is 10+
denominator <- sum(hh_sizes$prop)
hh_dist <- hh_sizes[, list(hh_size_prop=sum(prop)/denominator), by="hh_size"]
ten_plus <- sum(hh_dist[hh_size>=10]$hh_size_prop)
hh_dist <- hh_dist[hh_size<=10]
hh_dist[hh_size==10, hh_size_prop:=ten_plus]

if (sum(hh_dist$hh_size_prop)!=1){
  warning("Household size distribution improperly computed!")
}

# merge onto stock and flow values
stock_and_flow <- merge(stock_and_flow, hh_dist, by="hh_size", all=T)

# weight stock and flow values by household propotions 
stock_and_flow[, weighted_prob_no_nets:=hh_size_prop*SF_prob_no_nets]
stock_and_flow[, weighted_prob_any_net:=hh_size_prop*(1-SF_prob_no_nets)]

# calculate year-month-country access
stock_and_flow[, nat_access:=calc_access(hh_size, weighted_prob_no_nets, weighted_prob_any_net, SF_mean_nets_per_hh), 
                by=list(iso3, hh_size, year)]
stock_and_flow[, emplogit_nat_access:=emplogit(nat_access, 1000)] # todo: still don't understand this emplogit calc


## Create INLA Prediction objects  ## ---------------------------------------------------------

INLA:::inla.dynload.workaround() 

spatial_mesh = inla.mesh.2d(loc= unique(prediction_xyz[, list(x,y,z)]),
                                cutoff=0.006,
                                min.angle=c(25,25),
                                max.edge=c(0.06,500) )
temporal_mesh=inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2)

predictive_A <-inla.spde.make.A(spatial_mesh, 
                               loc=prediction_xyz[, list(x,y,z)], 
                               group=rep(min(this_year, 2015), length(prediction_indices)),
                               group.mesh=temporal_mesh ) 

