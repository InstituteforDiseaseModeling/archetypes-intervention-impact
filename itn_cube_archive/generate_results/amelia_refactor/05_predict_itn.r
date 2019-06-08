###############################################################################################################
## 05_predict_itn.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt's original code to predict access and use outputs from the previously-run
## input models. 
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor joint_dir=gs://map_data_z/users/amelia/itn_cube/input_data_archive/ --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/05_predict_itn.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam/05_predictions --command 'Rscript ${CODE}'

rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo", "VGAM", "raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam//"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam/05_predictions/"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data_archive"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
  cov_dir <- "/Volumes/GoogleDrive/Shared drives/cubes/5km incomplete/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
  cov_dir <- Sys.getenv("cov_dir")
}

# directories
set.seed(212)

# TODO: load shapefiles and templates only when you're actually plotting (objects named World, Africa, Water, img)
print("loading inla outputs and relevant functions")
# Load access and gap INLA outputs
load(file.path(input_dir, "03_access_deviation.Rdata"))
load(file.path(input_dir, "04_use_gap.Rdata"))

# load function script
source(file.path(func_dir, "01_create_database_functions.r")) # for calc_access
source(file.path(func_dir, "03_05_general_functions.r"))

# why does *this* go to 2016 when our inla only goes to 2014/15?
all_prediction_years <- 2000:2016

# temp
prediction_years <- 2000:2016

## Load Covariates  ## ---------------------------------------------------------


for (this_year in prediction_years){
  print(paste("predicting for year", this_year))
  
  print("loading covariates")
  static_covs <- fread(file.path(input_dir, "02_static_covariates.csv"))
  
  
  # TEMP TO COINCIDE WITH SAM: restrict the covariate year to 2001-2013
  cov_year <- min(2013, max(this_year, 2001))
  
  annual_covs <- fread(file.path(input_dir, "02_annual_covariates.csv"))
  annual_covs <- annual_covs[year==cov_year]
  dynamic_covs <- fread(file.path(input_dir, paste0("02_dynamic_covariates/dynamic_", cov_year, ".csv")))
  
  # more temp ---
  annual_covs[, year:= this_year]
  dynamic_covs[, year:=this_year]
  # ----
  
  annual_covs <- merge(static_covs, annual_covs, by="cellnumber", all=T)
  
  
  ## Get locations in x-y-z space of each pixel centroid for prediction ## ---------------------------------------------------------
  national_raster_dir <- file.path(joint_dir, "general/african_cn5km_2013_no_disputes.tif")
  national_raster <- raster(national_raster_dir)
  NAvalue(national_raster) <- -9999
  
  prediction_indices <- static_covs$cellnumber
  prediction_cells <- data.table(row_id=prediction_indices, gaul=extract(national_raster, prediction_indices))
  
  prediction_cells <- cbind(prediction_cells, data.table(xyFromCell(national_raster, prediction_indices)))
  setnames(prediction_cells, c("x", "y"), c("longitude", "latitude"))
  prediction_xyz <- ll_to_xyz(prediction_cells)
  
  
  
  ## Get national access means from stock and flow  ## ---------------------------------------------------------
  
  print("loading and calculating national access from stock and flow")
  # TODO: move all of this to prep_database, it doesn't need to be here
  orig_stock_and_flow <- fread(file.path(input_dir, "01_stock_and_flow_prepped.csv"))
  stock_and_flow <- orig_stock_and_flow[year>=this_year & year<(this_year+1)] # keep only the years we want to predict
  
  iso_gaul_map<-fread(file.path(joint_dir, "general/country_table_populations.csv")) # load table to match gaul codes to country names
  setnames(iso_gaul_map, c("GAUL_CODE", "COUNTRY_ID", "NAME"), c("gaul", "iso3", "country"))

  prediction_cells <- merge(prediction_cells, iso_gaul_map, by="gaul", all.x=T)
  
  # load household size distributions for each survey
  # TODO: rename this to "stock and flow", move national use vals to this folder
  hh_sizes<-fread(file.path(joint_dir, "stock_and_flow/HHsize.csv"))
  survey_key=fread(file.path(joint_dir, "stock_and_flow/hh_surveys_key.csv"))
  
  # format hh_sizes
  hh_sizes[, V1:=NULL]
  hh_sizes <- melt.data.table(hh_sizes, id.vars="HHSurvey", variable.name="hh_size", value.name="prop")
  hh_sizes[, hh_size:=as.integer(hh_size)]
  
  # BUG: The original script seems like it subsets on country-specific survey, 
  # but I think it always uses the entire dataset to calculate home distributions because it filters on ISO instead of full country name
  
  # update country names to fit with iso_gaul_map (for use when you actually DO calculate props by country name)
  survey_key[, Status:=NULL]
  setnames(survey_key, c("Svy Name", "Name"), c("HHSurvey", "country"))
  
  # add two surveys not present in key
  survey_key <- rbind(survey_key, data.table(HHSurvey=c("Kenya 2007", "Rwanda 2013"),
                                             country=c("Kenya", "Rwanda")))
  
  # merge with name map and survey distribution
  survey_key <- merge(survey_key, iso_gaul_map, by="country", all.x=T)
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
      predictions[, final_prediction := inv_ihs(full, theta=theta)] # TODO: confirm the inverse ihs function is correct
      
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
  
  print("finding and plotting access")
  setnames(acc_dev_predictions, "final_prediction", "access_deviation")
  acc_dev_predictions <- merge(acc_dev_predictions, stock_and_flow_access, by=c("iso3", "year", "month"), all.x=T)
  acc_dev_predictions[, emplogit_access:= emplogit_nat_access + access_deviation]
  
  # todo: definitely something weird with emplogit and access dev
  acc_dev_predictions_transformed <- acc_dev_predictions[, list(iso3, year, month, cellnumber, 
                                                    emplogit_access=emplogit_access,
                                                    access=plogis(emplogit_access),
                                                    access_deviation=plogis(access_deviation),
                                                    nat_access=plogis(emplogit_nat_access))]
  
  # TODO: should na.rm be T or F?
  summary_access <- acc_dev_predictions_transformed[, list(nat_access=mean(nat_access, na.rm=F),
                                               access_deviation=mean(access_deviation, na.rm=F),
                                               access=mean(access, na.rm=F)),
                                        by=list(iso3, year, cellnumber)]
  summary_access <- summary_access[order(cellnumber)]
  
  
  access_map <- copy(national_raster)
  access_map[summary_access$cellnumber] <- summary_access$access
  access_map[!is.na(national_raster) & is.na(access_map)] <- 0
  
  deviation_map <- copy(national_raster)
  deviation_map[summary_access$cellnumber] <- summary_access$access_deviation
  deviation_map[!is.na(national_raster) & is.na(deviation_map)] <- 0
  
  nat_mean_map <- copy(national_raster)
  nat_mean_map[summary_access$cellnumber] <- summary_access$nat_access
  nat_mean_map[!is.na(national_raster) & is.na(nat_mean_map)] <- 0
  
  # write files
  print("writing access tifs")
  writeRaster(access_map, file.path(output_dir, paste0("ITN_",this_year,".ACC.tif")),NAflag=-9999,overwrite=TRUE)
  writeRaster(deviation_map, file.path(output_dir, paste0("ITN_",this_year,".DEV.tif")),NAflag=-9999,overwrite=TRUE)
  writeRaster(nat_mean_map, file.path(output_dir, paste0("ITN_", this_year,".MEAN.tif")),NAflag=-9999,overwrite=TRUE)
  
  # DELETE AFTER MODEL COMPARE: write files with buggy transformations to mimic Sam's code for deviation and national mean
  sam_odd <- acc_dev_predictions[, list(nat_access=mean(emplogit_nat_access),
                                        access_deviation=mean(access_deviation)), by=list(iso3, year, cellnumber)]
  sam_odd <- sam_odd[, list(iso3, year, cellnumber, nat_access=plogis(nat_access),
                                access_deviation, 
                                access=summary_access$access)]
  
  odd_deviation_map <- copy(national_raster)
  odd_deviation_map[sam_odd$cellnumber] <- sam_odd$access_deviation
  odd_deviation_map[!is.na(national_raster) & is.na(odd_deviation_map)] <- 0
  
  odd_nat_mean_map <- copy(national_raster)
  odd_nat_mean_map[sam_odd$cellnumber] <- sam_odd$nat_access
  odd_nat_mean_map[!is.na(national_raster) & is.na(odd_nat_mean_map)] <- 0
  
  writeRaster(odd_deviation_map, file.path(output_dir, paste0("ITN_",this_year,".ODD_DEV.tif")),NAflag=-9999,overwrite=TRUE)
  writeRaster(odd_nat_mean_map, file.path(output_dir, paste0("ITN_", this_year,".ODD_MEAN.tif")),NAflag=-9999,overwrite=TRUE)
  
  ## Predict use  ## ---------------------------------------------------------
  
  
  print("predicting use gap")
  use_gap_predictions <- predict_inla(model=mod_pred_use, 
                                      theta=theta_use)
  
  print("finding and plotting use")
  setnames(use_gap_predictions, "final_prediction", "use_gap")
  
  use_gap_predictions <- merge(use_gap_predictions, acc_dev_predictions, by=c("iso3", "year", "month", "cellnumber"), all=T)
  use_gap_predictions[, emplogit_use:= emplogit_access - use_gap]
  
  use_gap_predictions_transformed <- use_gap_predictions[, list(iso3, year, month, cellnumber, 
                                                    use=plogis(emplogit_use),
                                                    use_gap=plogis(use_gap))]
  
  summary_use <- use_gap_predictions_transformed[, list(use=mean(use, na.rm=F),
                                            use_gap=mean(use_gap, na.rm=F)),
                                     by=list(iso3, year, cellnumber)]
  summary_use <- summary_use[order(cellnumber)]
  
  
  use_map <- copy(national_raster)
  use_map[summary_use$cellnumber] <- summary_use$use
  use_map[!is.na(national_raster) & is.na(use_map)] <- 0
  
  use_gap_map <- copy(national_raster)
  use_gap_map[summary_use$cellnumber] <- summary_use$use_gap
  use_gap_map[!is.na(national_raster) & is.na(use_gap_map)] <- 0
  
  # write files
  print("writing use tifs")
  writeRaster(use_map, file.path(output_dir, paste0("ITN_",this_year,".USE.tif")),NAflag=-9999,overwrite=TRUE)
  writeRaster(use_gap_map, file.path(output_dir, paste0("ITN_",this_year,".GAP.tif")),NAflag=-9999,overwrite=TRUE)
  
  
  # DELETE AFTER MODEL COMPARE: write files with buggy transformations to mimic Sam's code for deviation and national mean
  sam_odd_use <- use_gap_predictions[, list(use_gap=mean(use_gap)), by=list(iso3, year, cellnumber)]
  
  odd_gap_map <- copy(national_raster)
  odd_gap_map[sam_odd_use$cellnumber] <- sam_odd_use$use_gap
  odd_gap_map[!is.na(national_raster) & is.na(odd_gap_map)] <- 0
  
  writeRaster(odd_gap_map, file.path(output_dir, paste0("ITN_", this_year,".ODD_GAP.tif")),NAflag=-9999,overwrite=TRUE)
  
}


## "squash" certain years to zero  ## ---------------------------------------------------------

if (2016 %in% prediction_years){
  print("Squashing early years") # check: I think these are use estimates from the bucket model
  stock_and_flow_use <- fread(file.path(joint_dir, "stock_and_flow/indicators_access_qtr_new.csv"))
  stock_and_flow_use <- stock_and_flow_use[2:nrow(stock_and_flow_use)]
  names(stock_and_flow_use) <- c("country", seq(2000, 2018, 0.25))
  stock_and_flow_use <- melt(stock_and_flow_use, id.vars = "country", variable.name="year", value.name="use")
  stock_and_flow_use[, year:=floor(as.numeric(as.character(year)))]
  
  annual_use <- stock_and_flow_use[year<=max(all_prediction_years), list(use=mean(use)), by=list(country, year)]
  annual_use[, to_squash:=as.integer(use<0.02)]
  squash_map <- dcast.data.table(annual_use, country~year, value.var="to_squash")
  
  restrict_indicator <- function(prior_val, current_val){
    return(ifelse(prior_val==0, 0, current_val))
  }
  
  for (year in all_prediction_years[2:length(all_prediction_years)]){
    squash_map[[as.character(year)]] <- restrict_indicator(squash_map[[as.character(year-1)]], squash_map[[as.character(year)]])
  }
  
  squash_map <- melt(squash_map, id.vars = "country", value.name="to_squash", variable.name="year")
  squash_map[, year:=as.integer(as.character(year))]
  squash_map <- merge(squash_map, iso_gaul_map, by="country", all.x=T)
  
  
  # for all the country-years where "indicator" equals 1, set use to 0
  squash_map[, country_count:= sum(to_squash), by=list(year)] 
  
  years_to_squash <- unique(squash_map[country_count>0]$year)
  for (this_year in years_to_squash){
    
    # temp until you can figure out how to run everything in one batch
    if (file.exists(file.path(output_dir, paste0("ITN_", this_year, ".USE.tif")))){
      orig_use <- raster(file.path(output_dir, paste0("ITN_", this_year, ".USE.tif")))
    }else{
      orig_use <- raster(file.path(input_dir, "05_predictions", paste0("ITN_", this_year, ".USE.tif")))
    }
    
    new_use <- copy(orig_use)
    gauls_to_squash <- squash_map[year==this_year & to_squash==1]$gaul
    new_use[national_raster %in% gauls_to_squash] <- 0
    
    writeRaster(new_use, file.path(output_dir, paste0("ITN_",this_year,".RAKED_USE.tif")),NAflag=-9999,overwrite=TRUE)
  }
  
}


