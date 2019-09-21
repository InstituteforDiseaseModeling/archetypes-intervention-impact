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
out_dir <- file.path(base_dir, "no_transmission_limits")

# if synoptic==T, take the mean across years. Ohterwise save a raster for every year-month
synoptic <- T
overwrite_extraction <- T

template <- raster(file.path(in_dir, "MAP_Regions_Pf_5k.tif"))

input_files <- list.files(in_dir)[list.files(in_dir) %like% paste0("\\.nc$")]

# TEMP
input_files <- input_files[[1]]

for (this_input_file in input_files){
  
  cov_name <- gsub("\\.nc", "", this_input_file)
  print(paste("Extracting ERA5 for", cov_name))
  
  all_layers <- brick(file.path(in_dir, this_input_file))
  
  name_pattern <- "X([0-9]{4})\\.([0-9]{2}).*"
  names(all_layers) <- gsub(name_pattern, paste0(cov_name, "_year_", "\\1", "_month_", "\\2"), names(all_layers))
  
  if (synoptic){
    print("Finding synoptic values")
    synoptic_layers <- lapply(1:12, function(this_month){
      print(this_month)
      month_str <- str_pad(as.character(this_month), width=2, side="left", pad="0")
      stack_names <- names(all_layers)[names(all_layers) %like% paste0("month_", month_str)]
      
      this_stack <- all_layers[[stack_names]]
      synoptic_mean <- calc(this_stack, fun=mean)
      names(synoptic_mean) <- paste0(cov_name, "_month_", month_str)
      return(synoptic_mean)
    })
    final_layers <- stack(synoptic_layers)
  }else{
    final_layers <- all_layers
  }
  
  print("rotating and converting units")
  # rotate to -180, 180
  new_extent <- extent(c(0, 360, ymin(extent(final_layers)), ymax(extent(final_layers))))
  extent(final_layers) <- new_extent
  final_layers <- rotate(final_layers)
  
  if (cov_name %like% "air_temp"){
    # convert from kelvin to celsius
    final_layers <- final_layers - 273.15
  }else if (cov_name %like% "precip"){
    # convert from m to mm
    final_layers <- final_layers * 1000
  }
  
  print("resampling")
  # resample to match MAP resolution
  final_layers <- align_res(final_layers, template)
  
  # crop to each continent; save
  for (this_continent in c("africa", "americas", "asia")){
    print(paste("saving by continent", this_continent))
    this_out_dir <- file.path(out_dir, this_continent, cov_name)
    dir.create(this_out_dir, showWarnings=F, recursive=T)
    extraction_fname <- file.path(this_out_dir, paste0(cov_name, "_vals.csv"))
    
    if (file.exists(extraction_fname) & overwrite_extraction==F){
      print("values already extracted")
    }else{
      print("saving rasters")
      mask_layer <- raster(file.path(this_out_dir, "../mask.tif"))
      cropped_layers <- extract_values(final_layers, out_fname = file.path(this_out_dir, paste0(names(final_layers), ".tif")),
                                       mask = mask_layer)
      
      if(!synoptic){
        stop("TODO: THINK ABOUT HOW TO STRUCTURE A MULTIVARIABLE FILE")
      }else{
        vals <- as.matrix(cropped_layers)
        vals <- data.table(vals,
                           id=1:nrow(vals))
        names(vals) <- gsub(paste0(cov_name, "_month_"), "", names(vals)) 
        vals <- vals[!is.na(vals[["01"]])]
        vals <- melt(vals, id.vars = "id", variable.name="variable_val")
        vals <- vals[, list(cov=cov_name,
                            variable_name="month",
                            variable_val,
                            id,
                            value
        )]
        
        print("saving extracted values")
        write.csv(vals, file=extraction_fname, row.names=F)
      }
      
    }
    
  }
  
}

