###############################################################################################################
## 02_prep_covariates.r
## Amelia Bertozzi-Villa
## April 2019
## 
## A restructuring of Sam Bhatt's original code to prepare covariates for the ITN cube model.
## 
##############################################################################################################

# todo: change dsub so that output path is joint_data
# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/cubes_5km joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor --input database_fname=gs://map_data_z/users/amelia/itn_cube/create_database/output/ITN_final_clean_access_18March2019_COMPARE.csv CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/02_prep_covariates.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/joint_data --command 'Rscript ${CODE}'

rm(list=ls())

# package installation/loading
package_load <- function(package_list){
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/joint_data"
  # todo: make database output go to joint_data
  database_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/output/ITN_final_clean_access_18March2019_COMPARE.csv"
  input_dir <- "/Volumes/GoogleDrive/Team Drives/cubes/5km incomplete/"
  # input_dir <- "/Volumes/map_data/cubes/5km/"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/joint_data"
} else {
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  database_fname <- Sys.getenv("database_fname") # location of output file from generate_database_refactored.r
  input_dir <- Sys.getenv("input_dir") # here, location of covariate data 
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
  output_dir <- Sys.getenv("output_dir")
}


output_fname <- file.path(output_dir, 'covariates_27April2019.csv')

# TODO: function for extracting a raster stack and applying it to data
extract_values <- function(raster_list, dataset, names=""){
  
  cov_stack <- stack(raster_list)
  NAvalue(cov_stack)=-9999
  extracted_covs <- data.table(cov_stack[dataset$cellnumber])
  if (names!=""){
    names(extracted_covs) <- names
  }
  
  dataset <- cbind(dataset, extracted_covs)
  
  return(dataset)
}



# Load data from create_database.r, and list of covariates  ------------------------------------------------------------
data<-fread(database_fname)
data[, fulldate:=as.Date(as.yearmon(year))]

cov_details <- fread(file.path(func_dir, "covariate_key.csv"))
cov_details[, used_sam:= as.logical(used_sam)]
cov_details <- cov_details[used_sam==T]


### Static covariates: just extract and apply to all data  ----------------------------------------------------------------------------#######################  

static_fnames <- cov_details[type=="static", list(fname=file.path(input_dir, cov_name, fpath_append, fname))]

# TODO: creates 3 NA values. Explore.
data <- extract_values(static_fnames$fname, data)

### Annual covariates: extract and apply by year  ----------------------------------------------------------------------------#######################  

annual_covs <- cov_details[type=="year"]

# find base fnames-- allows for the possibility of a non-year looping variable such as landcover fraction
base_fnames <- lapply(annual_covs$cov_name, function(this_covname){
  this_cov <- annual_covs[cov_name==this_covname]
  if (this_cov$other_var==""){
    fnames <- data.table(cov_name=this_covname,
                         base_fname=file.path(input_dir, this_covname, this_cov$fpath_append, this_cov$fname),
                         year_start=this_cov$year_start,
                         year_end=this_cov$year_end,
                         colname=this_covname 
    )
  }else{ 
    fname_list <- lapply(this_cov$othervar_start:this_cov$othervar_end, function(newval){
      new_fname <- str_replace(this_cov$fname, this_cov$other_var, as.character(newval))
      return(file.path(input_dir, this_covname, this_cov$fpath_append, new_fname))
    })
    fnames <- data.table(cov_name=this_covname,
                         base_fname=unlist(fname_list),
                         year_start=this_cov$year_start,
                         year_end=this_cov$year_end,
                         colname=this_cov$othervar_start:this_cov$othervar_end)
    fnames[, colname:=paste0(cov_name, "_", colname)]
  }
  return(fnames)
})

base_fnames <- rbindlist(base_fnames)


# subset data years, extracting year-specific covariates for each
annual_list <- lapply(unique(data$flooryear), function(this_year){
  print(this_year)
  these_fnames <- copy(base_fnames)
  these_fnames[, year_to_use:=pmin(year_end, this_year)] # cap year by covariate availability
  these_fnames[, year_to_use:=pmax(year_to_use, year_start)] 
  these_fnames[, full_fname:=str_replace(base_fname, "YEAR", as.character(year_to_use))]
  
  subset <- extract_values(these_fnames$full_fname, data[flooryear==year], names=these_fnames$colname)
  return(subset)
})
data <- rbindlist(annual_list)


### Fully dynamic covariates: extract and apply by month and year  ----------------------------------------------------------------------------#######################  

dynamic_covs <- cov_details[type=="yearmon" & cov_name!="LST_night"]

dynamic_list <- lapply(unique(data$fulldate), function(this_date){
  this_month <- str_pad(month(this_date), 2, pad="0")
  this_year <- year(this_date)
  
  these_fnames <- copy(dynamic_covs)
  these_fnames[, year_to_use:=pmin(year_end, this_year)] # cap year by covariate availability
  these_fnames[, year_to_use:=pmax(year_to_use, year_start)] 
  these_fnames[, new_fname:=str_replace(fname, "YEAR", as.character(year_to_use))]
  these_fnames[, new_fname:=str_replace(new_fname, "MONTH", this_month)]
  these_fnames[, full_fname:=file.path(input_dir, cov_name, fpath_append, new_fname)]
  
  subset <- extract_values(these_fnames$full_fname, data[fulldate==this_date], names=these_fnames$cov_name) 
  
})
data <- rbindlist(dynamic_list)

write.csv(data, output_fname, row.names = F)

