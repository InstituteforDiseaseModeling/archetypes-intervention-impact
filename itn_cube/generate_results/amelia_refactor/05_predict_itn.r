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

package_load(c("zoo", "VGAM", "raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

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
print("loading inla outputs and relevant functions")
# Load access and gap INLA outputs
load(file.path(input_dir, "03_access_deviation.Rdata"))
load(file.path(input_dir, "04_use_gap.Rdata"))

# load function script
source(file.path(func_dir, "01_create_database_functions.r")) # for calc_access
source(file.path(func_dir, "03_access_dev_functions.r")) # for 1-val emplogit
source(file.path(func_dir, "05_predict_itn_functions.r"))

# why does *this* go to 2016 when our inla only goes to 2014/15?
prediction_years <- 2000:2016


## Load Covariates  ## ---------------------------------------------------------

this_year <- 2007

print(paste("predicting for year", this_year))

print("loading covariates")
static_covs <- fread(file.path(input_dir, "02_static_covariates.csv"))
annual_covs <- fread(file.path(input_dir, "02_annual_covariates.csv"))
annual_covs <- annual_covs[year==this_year]
dynamic_covs <- fread(file.path(input_dir, paste0("02_dynamic_covariates/dynamic_", this_year, ".csv")))

annual_covs <- merge(static_covs, annual_covs, by="cellnumber", all=T)


## Get locations in x-y-z space of each pixel centroid for prediction ## ---------------------------------------------------------
national_raster_dir <- file.path(joint_dir, "../african_cn5km_2013_no_disputes.tif")
national_raster <- raster(national_raster_dir)

prediction_indices <- static_covs$cellnumber
prediction_cells <- data.table(row_id=prediction_indices, gaul=extract(national_raster, prediction_indices))

prediction_cells <- cbind(prediction_cells, data.table(xyFromCell(national_raster, prediction_indices)))
setnames(prediction_cells, c("x", "y"), c("longitude", "latitude"))
prediction_xyz <- ll.to.xyz(prediction_cells)



## Get national access means from stock and flow  ## ---------------------------------------------------------

print("loading and calculating national access from stock and flow")
# TODO: move all of this to prep_database, it doesn't need to be here
orig_stock_and_flow <- fread(file.path(input_dir, "01_stock_and_flow_prepped.csv"))
stock_and_flow <- orig_stock_and_flow[year>=this_year & year<(this_year+1)] # keep only the years we want to predict
iso_gaul_map <- fread(file.path(joint_dir, "National_Config_Data.csv"))
iso_gaul_map <- iso_gaul_map[, list(iso3=ISO3, gaul=GAUL_Code, Name=IHME_Location_ASCII_Name)]

# add missing values 
iso_gaul_map <- rbind(iso_gaul_map, 
                      data.table(iso3="SHN",
                                 gaul=207,
                                 Name="Ascension and Tristan da Cunha")
                      )

# update cell file with iso3s (also fix odd gaul value for s africa)
prediction_cells[gaul==1013965, gaul:=227]
prediction_cells <- merge(prediction_cells, iso_gaul_map, by="gaul", all.x=T)

iso_gaul_map <- iso_gaul_map[iso3 %in% stock_and_flow$iso3]

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

stock_and_flow_access <- lapply(unique(stock_and_flow$iso3), function(this_iso){
  country_access <- lapply(unique(stock_and_flow$year), function(this_time){
    # print(paste(this_iso, ":", this_time))
    subset <- stock_and_flow[iso3==this_iso & year==this_time]
    access <- calc_access(subset, return_mean = T)
    return(data.table(iso3=this_iso, 
                      year=this_time,
                      nat_access=access)
                      )
  })
  return(rbindlist(country_access))
})

stock_and_flow_access <- rbindlist(stock_and_flow_access)
stock_and_flow_access[, emplogit_nat_access:=emplogit(nat_access, 1000)] # todo: still don't understand this emplogit calc
stock_and_flow_access <- merge(stock_and_flow_access, data.table(year=sort(unique(stock_and_flow_access$year)),
                                                                 month=1:12),
                               by="year", all=T)
stock_and_flow_access <- stock_and_flow_access[, list(iso3, year=this_year, month, nat_access, emplogit_nat_access)]


## Create INLA Prediction objects  ## ---------------------------------------------------------

print("creating INLA prediction objects")

INLA:::inla.dynload.workaround() 

# todo: why the capping at 2015?
spatial_mesh <-  copy(spatial_mesh_acc) # should be identical to spatial_mesh_use
temporal_mesh <- inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2)

A_matrix <-inla.spde.make.A(spatial_mesh, 
                               loc=as.matrix(prediction_xyz[, list(x,y,z)]), 
                               group=rep(min(this_year, 2015), length(prediction_indices)),
                               group.mesh=temporal_mesh ) 


## Predict access  ## ---------------------------------------------------------

predict_inla <- function(model, theta){
  fixed_effects <- model$summary.fixed
  random_effects <- model$summary.random$field
  predicted_random <- drop(A_matrix %*% random_effects$mean)
  
  all_predictions <- lapply(1:12, function(this_month){
    print(paste("predicting for month", this_month))
    these_covs <- merge(annual_covs, dynamic_covs[month==this_month], by=c("year", "cellnumber"), all=T)
    these_covs[, "Intercept":=1]
    
    predictions <- data.table(year=this_year,
                              month=this_month,
                              cellnumber=these_covs$cellnumber)
    
    predictions[, fixed:= as.matrix(these_covs[, rownames(fixed_effects), with=F]) %*% fixed_effects$mean]
    predictions[, random:= predicted_random]
    predictions[, full:= fixed + random]
    predictions[, access_deviation := Inv.IHS(full, theta=theta)] # TODO: confirm the inverse ihs function is correct
    
    predictions <- merge(predictions, prediction_cells[, list(cellnumber=row_id, iso3)], by="cellnumber", all=T)
    
    return(predictions) # for memory reasons
    # return(predictions[, list(year, month, cellnumber, final_prediction)])
  })
  
  all_predictions <- rbindlist(all_predictions)
  
  return(all_predictions)
}

print("predicting access deviation")
acc_dev_predictions <- predict_inla(model=mod_pred_acc, 
                                    theta=theta_acc)

acc_dev_predictions <- merge(acc_dev_predictions, stock_and_flow_access, by=c("iso3", "year", "month"), all.x=T)
acc_dev_predictions[, pixel_access:= emplogit_nat_access + access_deviation]


# todo: check what's going wrong in aggregation
summary_access <- acc_dev_predictions[, list(nat_access=plogis(mean(emplogit_nat_access)),
                                              access_deviation=plogis(mean(access_deviation)),
                                              access=plogis(mean(pixel_access))),
                                       by=list(iso3, year, cellnumber)]
summary_access <- summary_access[order(cellnumber)]

access_map <- copy(national_raster)
access_map[summary_access$cellnumber] <- summary_access$access
access_map[is.na(access_map)] <- 0
access_map[access_map==-9999] <- NA

deviation_map <- copy(national_raster)
deviation_map[summary_access$cellnumber] <- summary_access$deviation
deviation_map[is.na(deviation_map)] <- 0
deviation_map[deviation_map==-9999] <- NA




