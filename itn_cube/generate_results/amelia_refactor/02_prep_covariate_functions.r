###############################################################################################################
## 02_prep_covariates_functions.r
## Amelia Bertozzi-Villa
## May 2019
## 
## Functions to accompany 02_prep_covariates.r-- prepping survey data to go into the ITN cube model. 
## -which_non_null: way to find raster indices for which you want to extract other raster values 
## -extract_values: main function for loading a raster stack and returning the relevant values
## -get_annual_fnames: find filepaths for yearly covariates
##############################################################################################################

which_non_null <- function(raster_fname){
  reference_raster <- raster(raster_fname)
  NAvalue(reference_raster) <- -9999
  reference_vals <- getValues(reference_raster)
  index_vals <- which(!is.na(reference_vals))
  return(index_vals)
}

# function for extracting a raster stack and applying it to data
extract_values <- function(raster_fname_list, extraction_indices, names=c()){
  
  cov_stack <- stack(raster_fname_list)
  NAvalue(cov_stack)=-9999
  extracted_covs <- data.table(cov_stack[extraction_indices])
  if (length(names)>0){
    names(extracted_covs) <- names
  }
  
  return(extracted_covs)
}


get_annual_fnames <- function(covname, covariate_dt, input_dir){
  this_cov <- covariate_dt[cov_name==covname]
  if (this_cov$other_var==""){
    fnames <- data.table(cov_name=covname,
                         base_fname=file.path(input_dir, covname, this_cov$fpath_append, this_cov$fname),
                         start_year=this_cov$start_year,
                         end_year=this_cov$end_year,
                         colname=covname 
    )
  }else{ 
    fname_list <- lapply(this_cov$othervar_start:this_cov$othervar_end, function(newval){
      new_fname <- str_replace(this_cov$fname, this_cov$other_var, as.character(newval))
      return(file.path(input_dir, covname, this_cov$fpath_append, new_fname))
    })
    fnames <- data.table(cov_name=covname,
                         base_fname=unlist(fname_list),
                         start_year=this_cov$start_year,
                         end_year=this_cov$end_year,
                         colname=this_cov$othervar_start:this_cov$othervar_end)
    fnames[, colname:=paste0(cov_name, "_", colname)]
  }
  return(fnames)
}