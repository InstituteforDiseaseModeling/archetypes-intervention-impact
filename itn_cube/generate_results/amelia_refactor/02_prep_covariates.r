###############################################################################################################
## 02_prep_covariates.r
## Amelia Bertozzi-Villa
## April 2019
## 
## A restructuring of Sam Bhatt's original code to prepare covariates for the ITN cube model.
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/cubes_5km func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor --input database_fname=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor/01_database.csv CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/02_prep_covariates.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor/ --command 'Rscript ${CODE}'

rm(list=ls())

# package installation/loading
package_load <- function(package_list){
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  # todo: make database output go to joint_data
  database_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190503_total_refactor/01_database.csv"
  input_dir <- "/Volumes/GoogleDrive/Team Drives/cubes/5km incomplete/"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190503_total_refactor/"
} else {
  database_fname <- Sys.getenv("database_fname") # location of output file from generate_database_refactored.r
  input_dir <- Sys.getenv("input_dir") # here, location of covariate data 
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
  output_dir <- Sys.getenv("output_dir")
}


output_fname <- file.path(output_dir, "02_covariates.csv")

# function for extracting a raster stack and applying it to data
extract_values <- function(raster_list, dataset, names=c()){
  
  cov_stack <- stack(raster_list)
  NAvalue(cov_stack)=-9999
  extracted_covs <- data.table(cov_stack[dataset$cellnumber])
  if (length(names)>0){
    names(extracted_covs) <- names
  }
  
  dataset <- cbind(dataset, extracted_covs)
  
  return(dataset)
}


# Load data from create_database.r, and list of covariates  ------------------------------------------------------------
data<-fread(database_fname)
data[, fulldate:=as.Date(as.yearmon(year))]
data[, row_id:=row.names(data)]

cov_details <- fread(file.path(func_dir, "covariate_key.csv"))
cov_details[, used_sam:= as.logical(used_sam)]
cov_details <- cov_details[used_sam==T]


### Static covariates: just extract and apply to all data  ----------------------------------------------------------------------------#######################  

print("Extracting static covariates")
static_fnames <- cov_details[type=="static", list(fname=file.path(input_dir, cov_name, fpath_append, fname))]

# TODO: creates 3 NA values. Explore.
data <- extract_values(static_fnames$fname, data)
print("static covariates successfully extracted")

### Annual covariates: extract and apply by year  ----------------------------------------------------------------------------#######################  

print("Extracting annual covariates")
annual_covs <- cov_details[type=="year"]

# find base fnames-- allows for the possibility of a non-year looping variable such as landcover fraction
base_fnames <- lapply(annual_covs$cov_name, function(this_covname){
  this_cov <- annual_covs[cov_name==this_covname]
  if (this_cov$other_var==""){
    fnames <- data.table(cov_name=this_covname,
                         base_fname=file.path(input_dir, this_covname, this_cov$fpath_append, this_cov$fname),
                         start_year=this_cov$start_year,
                         end_year=this_cov$end_year,
                         colname=this_covname 
    )
  }else{ 
    fname_list <- lapply(this_cov$othervar_start:this_cov$othervar_end, function(newval){
      new_fname <- str_replace(this_cov$fname, this_cov$other_var, as.character(newval))
      return(file.path(input_dir, this_covname, this_cov$fpath_append, new_fname))
    })
    fnames <- data.table(cov_name=this_covname,
                         base_fname=unlist(fname_list),
                         start_year=this_cov$start_year,
                         end_year=this_cov$end_year,
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
  these_fnames[, year_to_use:=pmin(end_year, this_year)] # cap year by covariate availability
  these_fnames[, year_to_use:=pmax(year_to_use, start_year)] 
  these_fnames[, full_fname:=str_replace(base_fname, "YEAR", as.character(year_to_use))]
  
  subset <- extract_values(these_fnames$full_fname, data[flooryear==this_year], names=these_fnames$colname)
  return(subset)
})
data <- rbindlist(annual_list)

# remove Urban landcover for high collinearity with population
data[, IGBP_Landcover_13:=NULL]


print("annual covariates successfully extracted")

### Fully dynamic covariates: extract and apply by month and year  ----------------------------------------------------------------------------#######################  

print("Extracting dynamic covariates")

dynamic_covs <- cov_details[type=="yearmon"]

dynamic_list <- lapply(unique(data$fulldate), function(this_date){
  
  print(this_date)
  
  this_month <- month(this_date)
  this_year <- year(this_date)
  
  these_fnames <- copy(dynamic_covs)
  these_fnames[, year_to_use:=pmin(end_year, this_year)] # cap year by covariate availability
  these_fnames[, year_to_use:=pmax(year_to_use, start_year)] 
  # cap to make sure the specific month_year is available (2000 and 2014 don't have all months)
  these_fnames[end_year==year_to_use & end_month<this_month, year_to_use:=end_year-1]
  these_fnames[start_year==year_to_use & start_month>this_month, year_to_use:=start_year+1]
  
  these_fnames[, new_fname:=str_replace(fname, "YEAR", as.character(year_to_use))]
  these_fnames[, new_fname:=str_replace(new_fname, "MONTH", str_pad(this_month, 2, pad="0"))]
  these_fnames[, full_fname:=file.path(input_dir, cov_name, fpath_append, new_fname)]
  
  subset <- extract_values(these_fnames$full_fname, data[fulldate==this_date], names=these_fnames$cov_name) 
  
})
data <- rbindlist(dynamic_list)
data <- data[order(row_id)]

print("dynamic covariates successfully extracted")

print("saving to file")
write.csv(data, output_fname, row.names = F)


