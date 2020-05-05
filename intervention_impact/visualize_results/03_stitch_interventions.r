###############################################################################################################
## 03_stitch_interventions.r
## Amelia Bertozzi-Villa
## April 2020
## 
##############################################################################################################

library(data.table)
library(raster)
library(rasterVis)

rm(list=ls())

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results")
setwd(func_dir)
# source("pr_to_r0.r")
source("map_ii_functions.r")

analysis_subdir <- "20200426_int_history"
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
cluster_raster_dir <- file.path(base_dir, "../archetypes/results")
main_dir <- file.path(base_dir, analysis_subdir)
suffix <- "recur_out"
out_dir <- file.path(main_dir,"results", "rasters")
dir.create(out_dir, recursive = T, showWarnings = F)

africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"
raster_input_dir <- "~/Desktop/covid_mitigation_project/pfpr_mean_rasters_20200214"

africa_shp <- readOGR(africa_shp_dir)
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)


intervention_category_map <- raster(file.path(main_dir, "input", "intervention_history_categories.tif"))
all_intervention_maps <- stack(file.path(out_dir, paste0("pfpr_Africa",suffix, ".tif")))
comparison_map <- raster(file.path(raster_input_dir, "pr_2019_rmean_Africa.tif"))

all_intervention_maps <- all_intervention_maps[[2:nlayers(all_intervention_maps)]] # the first layer is the "baseline" surface from which the others are drawn
all_intervention_maps <- crop(all_intervention_maps, intervention_category_map)
all_intervention_maps <- raster::mask(all_intervention_maps, intervention_category_map)

stitched_map <- lapply(1:nlayers(all_intervention_maps), function(layer_idx){
  print(layer_idx)
  this_layer <- all_intervention_maps[[layer_idx]]
  this_mask <- intervention_category_map==layer_idx
  pr_masked <- raster::mask(this_layer, this_mask, maskvalue=0)
  return(pr_masked)
})

stitched_map <- do.call(merge, stitched_map)
names(stitched_map) <- "EMOD PfPR"

comparison_map <- raster::mask(comparison_map, stitched_map)
names(comparison_map) <- "MAP PfPR"

diff <- comparison_map-stitched_map


levelplot(diff, par.settings=rasterTheme(rev(brewer.pal(11, "PuOr"))), at=seq(-1, 1, 0.02),
          xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), 
          main="Difference between MAP and Intervention Impact, 2019",
          maxpixesl=2e5) + 
  latticeExtra::layer(sp.polygons(africa_shp))


