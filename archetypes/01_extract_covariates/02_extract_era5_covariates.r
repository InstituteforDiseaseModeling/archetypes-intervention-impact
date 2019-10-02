## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 02_extract_era5_covariates.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# September 2019
# 
# For temperature and rainfall, this script extracts covariate-specific
# data from ERA5 netCDF files, transforms them into a panel dataset, and saves it for further clustering analyses.  
# 
# For a detailed project write-up see
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------

library(stringr)
library(data.table)
library(rasterVis)
library(raster)
library(ncdf4)

rm(list=ls())

source("extraction_functions.r")

root_dir <- Sys.getenv("HOME")
base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/covariates")
in_dir <- file.path(base_dir, "era5_raw_data")
template <- raster(file.path(in_dir, "MAP_Regions_Pf_5k.tif"))
input_files <- list.files(in_dir)[list.files(in_dir) %like% paste0("\\.nc$")]

### INPUTS: change these as desired #####----------------------------------------------------------------------------------------------------------------------------------

# please make your own output folder to avoid overwriting someone's work
out_dir <- file.path(base_dir, "example_for_jaline")

# set years_to_keep to NULL if you would like to keep all years
years_to_keep <- 2014:2015

# if synoptic==T, take the mean across years. Otherwise save a raster for every year-month
synoptic <- F

# if resample==T, change the resolution and extent of the rasters to match the "template" object above
resample <- F

# if you add any continent other than "global", you must have a raster named "mask.tif" in each 
# of your output directories. The global raster is cropped and masked to the "mask" file. 
continents_to_keep <- c("global")

# NOT RECOMMENDED IF SYNOPTIC==F: writes a big file
# if save_to_csv==T, all the values in the rasters will also be saved as a .csv in long format, 
# one per input file, where "id" refers to the cell number to which that value corresponds
save_to_csv <- F

### Main extraction loop #####----------------------------------------------------------------------------------------------------------------------------------

for (this_input_file in input_files){
  
  # load raw ERA5 data
  cov_name <- gsub("\\.nc", "", this_input_file)
  print(paste("Extracting ERA5 for", cov_name))
  
  all_layers <- brick(file.path(in_dir, this_input_file))
  
  name_pattern <- "X([0-9]{4})\\.([0-9]{2}).*"
  names(all_layers) <- gsub(name_pattern, paste0(cov_name, "_year_", "\\1", "_month_", "\\2"), names(all_layers))
  
  # subset to the years of interest, if supplied
  if (length(years_to_keep)>0){
    names_to_keep <- unlist(lapply(years_to_keep, function(year){
      names(all_layers)[names(all_layers) %like% year]
    }))
    all_layers <- all_layers[[names_to_keep]]
  }
  
  # take the mean across years, if requested
  if (synoptic){
    print("Finding synoptic values")
    synoptic_layers <- lapply(1:12, function(this_month){
      print(this_month)
      month_str <- str_pad(as.character(this_month), width=2, side="left", pad="0")
      stack_names <- names(all_layers)[names(all_layers) %like% paste0("month_", month_str)]
      
      this_stack <- all_layers[[stack_names]]
      synoptic_mean <- max(calc(this_stack, fun=mean), 0) # find mean, remove negative values
      names(synoptic_mean) <- paste0(cov_name, "_month_", month_str)
      return(synoptic_mean)
    })
    final_layers <- stack(synoptic_layers)
  }else{
    final_layers <- all_layers
  }
  
  # rotate from 0,360 to -180, 180
  print("rotating")
  new_extent <- extent(c(0, 360, ymin(extent(final_layers)), ymax(extent(final_layers))))
  extent(final_layers) <- new_extent
  final_layers <- rotate(final_layers)
  
  print("converting units")
  if (cov_name %like% "air_temp"){
    # convert from kelvin to celsius
    final_layers <- final_layers - 273.15
  }else if (cov_name %like% "precip"){
    # each monthly value should be multiplied by the number of days in the month
    month_durations <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    
    # extend to the length of all_layers. NOTE: this won't work correctly if you year doesn't start in January
    if (!names(final_layers)[[1]] %like% "month_01"){
      stop("Your time series doesn't start in January! Adjust the values of month_durations")
    }
    month_durations <- rep(month_durations, round(nlayers(final_layers)/12))[1:nlayers(final_layers)]
    
    # convert from m/day to mm/month
    final_layers <- final_layers * month_durations * 1000
  }
  
  # move to "template" resolution, if desired 
  if (resample){
    print("resampling")
    # resample to match MAP resolution
    final_layers <- align_res(final_layers, template)
  }
  
  
  # crop to each continent; save
  for (this_continent in continents_to_keep){
    print(paste("saving for continent:", this_continent))
    this_out_dir <- file.path(out_dir, this_continent, cov_name)
    dir.create(this_out_dir, showWarnings=F, recursive=T)
    extraction_fname <- file.path(this_out_dir, paste0(cov_name, "_vals.csv"))
    
    # save rasters
    if (this_continent=="global"){
      writeRaster(final_layers, file.path(this_out_dir, paste0(names(final_layers), ".tif")), bylayer=T, overwrite=T)
      vals <- as.matrix(final_layers)
    }else{
      mask_layer <- raster(file.path(this_out_dir, "../mask.tif"))
      cropped_layers <- crop_raster(final_layers, mask = mask_layer, 
                                    out_fname = file.path(this_out_dir, paste0(names(final_layers), ".tif"))
                                    )
      vals <- as.matrix(cropped_layers)
      rm(cropped_layers); gc()
    }
    
    # save as .csv, if desired
    if(save_to_csv){
      vals <- data.table(vals,
                         id=1:nrow(vals))
      print("reshaping and formatting columns")
      if(!synoptic){
        vals <- melt(vals, id.vars = "id", variable.name="variable_val")
        vals <- vals[, list(cov=cov_name,
                            year=as.integer(gsub(".*_year_([0-9]{4}).*", "\\1", variable_val)),
                            month=as.integer(gsub(".*_month_([0-9]{2})", "\\1", variable_val)),
                            id,
                            value
        )]
      }else{
        names(vals) <- gsub(paste0(cov_name, "_month_"), "", names(vals)) 
        vals <- vals[!is.na(vals[["01"]])]
        vals <- melt(vals, id.vars = "id", variable.name="variable_val")
        vals <- vals[, list(cov=cov_name,
                            variable_name="month",
                            variable_val,
                            id,
                            value
        )]
        
      }
      print("saving extracted values")
      write.csv(vals, file=extraction_fname, row.names=F)
    }
    
  } 
  rm("vals","all_layers", "final_layers"); gc()
}

