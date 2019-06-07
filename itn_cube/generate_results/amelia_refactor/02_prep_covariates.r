###############################################################################################################
## 02_prep_covariates.r
## Amelia Bertozzi-Villa
## April 2019
## 
## A restructuring of Sam Bhatt's original code to prepare covariates for the ITN cube model.
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/cubes_5km func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor joint_dir=gs://map_data_z/users/amelia/itn_cube/input_data_archive/ --input database_fname=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam/01_database.csv CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/02_prep_covariates.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam/ --command 'Rscript ${CODE}'

rm(list=ls())

# package installation/loading
package_load <- function(package_list){
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  database_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam//01_database.csv"
  input_dir <- "/Volumes/GoogleDrive/Shared drives/cubes/5km incomplete/"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data_archive"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam//"
} else {
  database_fname <- Sys.getenv("database_fname") # location of output file from generate_database_refactored.r
  input_dir <- Sys.getenv("input_dir") # here, location of covariate data 
  joint_dir <- Sys.getenv("joint_dir")
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
  output_dir <- Sys.getenv("output_dir")
}

prediction_years <- 2000:2016

source(file.path(func_dir, "02_prep_covariate_functions.r"))

# Load data from create_database.r, and list of covariates  ------------------------------------------------------------
data<-fread(database_fname)
data[, fulldate:=as.Date(as.yearmon(year))]
data[, month:=month(fulldate)]
data[, row_id:=as.integer(row.names(data))]

cov_dt <- fread(file.path(func_dir, "covariate_key.csv"))
cov_dt[, used_sam:= as.logical(used_sam)]
cov_dt <- cov_dt[used_sam==T]

# find the "valid" cell values for which we want to predict in step 5
raster_indices <- which_non_null(file.path(joint_dir, "general/african_cn5km_2013_no_disputes.tif"))


### Static covariates  ----------------------------------------------------------------------------#######################  

print("Extracting static covariates")
static_fnames <- cov_dt[type=="static", list(fname=file.path(input_dir, cov_name, fpath_append, fname))]

all_static <- extract_values(static_fnames$fname,raster_indices)
all_static[, cellnumber:=raster_indices]
write.csv(all_static, file.path(output_dir, "02_static_covariates.csv"), row.names = F)

data <- merge(data, all_static, by="cellnumber", all.x=T)

print("static covariates successfully extracted")

### Annual covariates  ----------------------------------------------------------------------------#######################  

print("Extracting annual covariates")
annual_cov_dt <- cov_dt[type=="year"]

# find base fnames-- allows for the possibility of a non-year looping variable such as landcover fraction
base_fnames <- lapply(annual_cov_dt$cov_name, get_annual_fnames, covariate_dt=annual_cov_dt, input_dir=input_dir)
base_fnames <- rbindlist(base_fnames)
base_fnames <- base_fnames[colname!="IGBP_Landcover_13"] # remove Urban landcover for high collinearity with population

print("Extracting whole-continent values")
ncores <- detectCores()
print(paste("--> Machine has", ncores, "cores available"))
registerDoParallel(ncores-2)

all_annual <- foreach(this_year=prediction_years) %dopar%{
  
  print(this_year)
  
  these_fnames <- copy(base_fnames)
  these_fnames[, year_to_use:=pmin(end_year, this_year)] # cap year by covariate availability
  these_fnames[, year_to_use:=pmax(year_to_use, start_year)] 
  these_fnames[, full_fname:=str_replace(base_fname, "YEAR", as.character(year_to_use))]
  
  subset <- extract_values(these_fnames$full_fname, raster_indices, names=these_fnames$colname)
  subset[, year:=this_year]
  subset[, cellnumber:=raster_indices]
  setcolorder(subset, c("year", "cellnumber", these_fnames$colname))
  
  # TEMP: name to correlate to old dataset
  new_names <- str_replace(names(subset), "IGBP_Landcover_", "landcover")
  new_names <- str_replace(new_names, "AfriPop", "populatopn")
  names(subset) <- new_names
  
  return(subset)
}

all_annual <- rbindlist(all_annual)


# isolate values for data
setnames(all_annual, "year", "flooryear")
data <- merge(data, all_annual, by=c("flooryear", "cellnumber"), all.x=T)
setnames(all_annual, "flooryear", "year")
write.csv(all_annual, file.path(output_dir, "02_annual_covariates.csv"), row.names = F)

print("annual covariates successfully extracted")


### Fully dynamic covariates: extract and apply by month and year  ----------------------------------------------------------------------------#######################  

print("Extracting dynamic covariates")

dynamic_cov_dt <- cov_dt[type=="yearmon"]
all_yearmons <- data.table(expand.grid(1:12, prediction_years))
names(all_yearmons) <- c("month", "year")

registerDoParallel(ncores-2)

dynamic_outdir <- file.path(output_dir, "02_dynamic_covariates")
if (!dir.exists(dynamic_outdir)){
  dir.create(dynamic_outdir)
}

all_dynamic <- foreach(month_index=1:nrow(all_yearmons)) %dopar% {
  
  print(all_yearmons[month_index])
  
  this_month <- all_yearmons[month_index]$month
  this_year <- all_yearmons[month_index]$year
  
  these_fnames <- copy(dynamic_cov_dt)
  these_fnames[, year_to_use:=pmin(end_year, this_year)] # cap year by covariate availability
  these_fnames[, year_to_use:=pmax(year_to_use, start_year)] 
  # cap to make sure the specific month_year is available (2000 and 2014 don't have all months)
  these_fnames[end_year==year_to_use & end_month<this_month, year_to_use:=end_year-1]
  these_fnames[start_year==year_to_use & start_month>this_month, year_to_use:=start_year+1]
  
  these_fnames[, new_fname:=str_replace(fname, "YEAR", as.character(year_to_use))]
  these_fnames[, new_fname:=str_replace(new_fname, "MONTH", str_pad(this_month, 2, pad="0"))]
  these_fnames[, full_fname:=file.path(input_dir, cov_name, fpath_append, new_fname)]
  
  subset <- extract_values(these_fnames$full_fname, raster_indices, names=these_fnames$cov_name)
  subset[, year:=this_year]
  subset[, month:=this_month]
  subset[, cellnumber:=raster_indices]
  setcolorder(subset, c("year", "month", "cellnumber", these_fnames$cov_name))
  
  # TEMP: name to correlate to old dataset
  setnames(subset, c("EVI", "LST_day", "LST_night", "TCW", "TSI"), c("evy", "lst_day", "lst_night", "tcw", "tsi"))
  ###
  
  return(subset)
}
all_dynamic <- rbindlist(all_dynamic)

# isolate values for data
setnames(all_dynamic, "year", "flooryear")


# -------TEMP: for the four cell values that are mis-extracted, preserve the bug for now by replacing the correct with the incorrect values 
incorrect_points <- data.table(cellnumber=c(901138, 1027886, 1603132, 2141192),
                               flooryear=c(2008, 2014, 2014, 2004),
                               month=c(12, 2, 10, 12),
                               sub_cellnumber=c(896098, 1026206, 1591372, 2137832))

data <- merge(data, incorrect_points, by=c("cellnumber", "flooryear", "month"), all.x=T)
data[is.na(sub_cellnumber), sub_cellnumber:=cellnumber]

setnames(all_dynamic, "cellnumber", "sub_cellnumber")
data <- merge(data, all_dynamic, by=c("flooryear", "month", "sub_cellnumber"), all.x=T)
data[, sub_cellnumber:=NULL]
setnames(all_dynamic, "sub_cellnumber", "cellnumber")

### -----------------------------------

setnames(all_dynamic, "flooryear", "year")

# save dynamic covariates(by year)
print("saving dynamic covariates by year")
for (this_year in prediction_years){
  write.csv(all_dynamic[year==this_year], file.path(dynamic_outdir, paste0("dynamic_", this_year,".csv")), row.names=F)
}

# write.csv(all_dynamic, file.path(output_dir, "02_dynamic_covariates.csv"), row.names=F)

# reorder
data <- data[order(row_id)]

print("dynamic covariates successfully extracted")

print("saving to file")
write.csv(data, file.path(output_dir, "02_data_covariates.csv"), row.names = F)





