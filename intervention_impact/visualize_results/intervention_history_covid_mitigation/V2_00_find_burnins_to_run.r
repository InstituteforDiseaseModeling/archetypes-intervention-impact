###############################################################################################################
## V2_00_find_burnins_to_run.r
## Amelia Bertozzi-Villa
## April 2020
## 
## Using a given set of burnins and initial raster, select only the relevant burnins to run (omit those that 
## are too high/low)
##############################################################################################################

library(raster)
library(data.table)

analysis_subdir <- "20200426_int_history"
burnin_subdir <- "20191009_mega_era5_new_arch"
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
out_dir <- file.path(base_dir, analysis_subdir, "input")

pfpr_dir <- "~/Desktop/covid_mitigation_project/pfpr_mean_rasters_20200214/"


# read burnin prevalences
burnin_vals <- fread(file.path(base_dir, burnin_subdir, "results/raw/", paste0("MAP_", burnin_subdir, "_Burnin.csv")))
burnin_vals <- burnin_vals[Site_Name>1 & Site_Name<11] # not interested in Asia or Americas
# max/mins across sites
summary_burnin_vals <- burnin_vals[, list(min_prev=min(initial_prev),
                                          max_prev=max(initial_prev),
                                          mean_prev=mean(initial_prev)),
                                   by=list(x_Temporary_Larval_Habitat)]
summary_burnin_vals <- summary_burnin_vals[order(x_Temporary_Larval_Habitat)]





# read and mask rasters
int_history_categories <- raster(file.path(out_dir, "intervention_history_categories.tif"))
no_int_prev <- raster(file.path(pfpr_dir, "pr_2000_rmean_Africa.tif"))
int_prev <- raster(file.path(pfpr_dir, "pr_2016_rmean_Africa.tif"))

no_int_prev <- raster::mask(no_int_prev, int_history_categories)
int_prev <- raster::mask(int_prev, int_history_categories)







