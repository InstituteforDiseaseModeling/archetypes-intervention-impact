## -----------------------------------------------------------------------------------------------------------------
# Quick raster comparison
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# September 2018
# 
# Ensure that previously saved vector abundance rasters are equivalent to the version I'd like to move to
## -----------------------------------------------------------------------------------------------------------------------

library(data.table)
library(raster)

rm(list=ls())

root_dir_local <- Sys.getenv("HOME")
root_dir_mount <- "/Volumes/map_data/mastergrids"

main_dir <- file.path(root_dir_local, "repos/malaria-atlas-project/seasonal_classification") 
setwd(main_dir)

covs <- fread("clustering_covariates.csv")
covs <- covs[variable=="species"]

species_names <- c("arabiensis", "funestus", "gambiae")

old_rasters <- lapply(species_names, function(species_name){
  this_cov <- covs[cov=="vector_abundance_old"]
  this_raster <- raster(file.path(root_dir_local, this_cov$dir, paste0(species_name, ".tif")))
  return(this_raster)
})
names(old_rasters) <- species_names

new_rasters <- lapply(c("arabiensis", "funestus_subgroup", "coluzzii_gambiae"), function(species_name){
  this_cov <- covs[cov=="relative_vector_abundance"]
  this_raster <- raster(file.path(root_dir_mount, this_cov$dir, paste0("Anopheles_", species_name, ".RelativeAbundance.tif")))
  return(this_raster)
})
names(new_rasters) <- species_names


diffs <- lapply(species_names, function(species_name){
  diff <- abs(old_rasters[[species_name]] - new_rasters[[species_name]])
  return(diff)
})





