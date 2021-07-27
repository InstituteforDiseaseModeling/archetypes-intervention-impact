###############################################################################################################
## 02_map_intervention_impact.r
## Amelia Bertozzi-Villa
## December 2019
## 
##############################################################################################################

library(data.table)
library(raster)
library(rasterVis)

rm(list=ls())

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results")
setwd(func_dir)
source("pr_to_r0.r")
source("map_ii_functions.r")

analysis_subdir <- "20191008_replicate_megatrends"
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
cluster_raster_dir <- file.path(base_dir, "../archetypes/results")
main_dir <- file.path(base_dir, analysis_subdir)
suffix <- ""
out_dir <- file.path(main_dir,"results", "rasters")
dir.create(out_dir, recursive = T, showWarnings = F)

africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"
raster_input_dir <- file.path(base_dir, "map_rasters/megatrends_project")

### Input variables  -----------------------------------------------------
plot_results <- F

calculate_par <- T
pop_fname <- file.path(raster_input_dir, "ssp5_total_2050_MG_5K.tif")

calculate_repro_number <- F

region <- "Africa" # geography for which you want to predict. currently responds only to "Africa" 

# where is your lookup table? 
lut_fname <- file.path(main_dir, paste0("results/clean/summary_impact", suffix, ".csv"))

# what interventions from the lookup table do you want to run? Select "all" if you want every intervention in your LUT to run
interventions <- "all"

# what rasters do you want to use as baselines? You'll get a new set of results for each
baseline_raster_fnames <- list( # "True PfPR 2017"="PfPR_rmean_Global_admin0_2017.tif",
 #  "MAP PfPR 2000"="pr_2000_rmean_Africa.tif"
  "Megatrends Base 2016"="actual_ssp2_base2000_2050.tif"
)

# What raster do you want to use as a maximum value on all results?
bounding_fname <- file.path(raster_input_dir, "MODEL43.2015.PR.ALL.rmean.tif")

# Are there additional rasters you want to visualize, but not apply a lookup table to?
comparison_fnames <- c()

### Shapefile  -----------------------------------------------------
africa_shp <- readOGR(africa_shp_dir)
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)

### Desired Geography  -----------------------------------------------------

site_details <- fread(file.path(main_dir, "input", "site_details.csv"))
if (region!="Global"){
  site_details <- site_details[continent==region]
}

# read in cluster rasters, assign values from LUT
cluster_map <- raster(file.path(cluster_raster_dir, unique(site_details$source), 
                                tolower(region), "02_kmeans",
                                paste0("map_", nrow(site_details),"_cluster.tif")))

### Lookup table -----------------------------------------------------

# read in lookup table
lut <- fread(lut_fname)

# keep only the desired interventions; simplify
if (!interventions=="all"){
  lut <- lut[int_id %in% interventions]
}else{
  interventions <- unique(lut$int_id)
}
lut <- unique(lut[, list(Site_Name, int_id, label, mean_initial, mean_final)])

# merge on a "cluster id" corresponding to the cluster values in each raster.
# For procedurally generated sites, this should be identical to Site_Name
# todo: update this with new cluster-id naming

# this will also drop the sites that are outside of the region of interest.
lut <- merge(lut, site_details[, list(Site_Name=id, cluster_id=cluster)], by="Site_Name")


### Align cluster and int history rasters -----------------------------------------------------
bounding_pr <- raster(bounding_fname)
cluster_map <- crop(cluster_map, bounding_pr)
cluster_map <- extend(cluster_map, bounding_pr)
cluster_map <- raster::mask(cluster_map, bounding_pr)

bounding_pr <- raster::mask(bounding_pr, cluster_map)

### Format any preliminary rasters  #####----------------------------------------------------------------------------------------------------------------------------------

# initialize list
pr_list <- list()
pr_idx <- 1

if (length(comparison_fnames)>0){
  print("Appending additional rasters")
  for (label in names(comparison_fnames)){
    print(label)
    this_raster <- raster(comparison_fnames[[label]])
    this_raster <- extend(crop(this_raster, bounding_pr), bounding_pr)
    
    raster_list[[raster_idx]] <- this_raster
    raster_idx <- raster_idx + 1
  }
}

### Main loop: apply lookup table  #####----------------------------------------------------------------------------------------------------------------------------------

for (baseline_label in names(baseline_raster_fnames)){
  print(baseline_label)
  baseline_pr <- raster(file.path(raster_input_dir, baseline_raster_fnames[[baseline_label]]))
  baseline_pr <- crop(baseline_pr, bounding_pr)
  baseline_pr <- extend(baseline_pr, bounding_pr)
  baseline_pr <- raster::mask(baseline_pr, cluster_map)
  names(baseline_pr) <- baseline_label
  
  pr_list[[pr_idx]] <- baseline_pr
  pr_idx <- pr_idx + 1 
  
  for (intervention in interventions){
    print(intervention)
    pr_list[[pr_idx]] <- apply_lookup(intervention, baseline_pr, lut, cluster_map, bounding_pr, bound=F)
    pr_idx <- pr_idx + 1 
  }
}

stacked_pr <- stack(pr_list)
rm(pr_list); gc();

# maybe: calculate all interventions as usual, and then stitch them together later!!!

# if (region=="Africa"){
#   stacked_pr <- extend(stacked_pr, africa_shp)
# }

# Plot and save
color_vals <- generate_full_pal()
int_count <- length(interventions)

if (plot_results){
  pdf(file.path(out_dir, paste0("pfpr_", region, ".pdf")), width=12, height=10)
  
  for(base_val in 1:length(baseline_raster_fnames)){
    end_idx <- (int_count+1) * base_val
    start_idx <- end_idx - int_count
    print(levelplot(stacked_pr[[start_idx:end_idx]], par.settings=rasterTheme(color_vals$pal), at=color_vals$breaks,
                    xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F)) +
            latticeExtra::layer(sp.polygons(africa_shp))
    )
  }
  
  graphics.off()
}

writeRaster(stacked_pr, options="INTERLEAVE=BAND", bylayer=F, suffix="names", filename=file.path(out_dir, paste0("pfpr_",region, suffix, ".tif")), overwrite=T)


### Convert to populations-at-risk  #####----------------------------------------------------------------------------------------------------------------------------------
if (calculate_par==T){
  
  print("Converting to Population-at-risk")
  
  pop <- raster(pop_fname)
  pop <- crop(pop, stacked_pr)
  
  has_transmission <- lapply(1:nlayers(stacked_pr), function(idx){
    this_layer <- copy(stacked_pr[[idx]])
    this_layer[this_layer<0.01] <- 0
    this_layer[this_layer>0] <- 1
    return(this_layer)
  })
  has_transmission <- stack(has_transmission)
  
  par <- has_transmission * pop
  names(par) <-  names(stacked_pr)
  par_millions <- cellStats(par, sum)/1000000
  
  # Plot and save
  if (plot_results)
  {
    pal <- wpal("seaside")
    pdf(file.path(out_dir, paste0("par_", region, ".pdf")), width=12, height=10)
    
    for(base_val in 1:length(baseline_raster_fnames)){
      end_idx <- (int_count+1) * base_val
      start_idx <- end_idx - int_count
      print(levelplot(par[[start_idx:end_idx]] + 1, par.settings=rasterTheme(pal), zscaleLog=T, margin=F, 
                      xlab=NULL, ylab=NULL, scales=list(draw=F))  +
              latticeExtra::layer(sp.polygons(africa_shp))
      )
    }
    
    graphics.off()
    
  }
  
  writeRaster(par, options="INTERLEAVE=BAND", bylayer=F, suffix="names", filename=file.path(out_dir, paste0("par_",region, suffix, ".tif")), overwrite=T)
  
}

### Convert to reproductive number  #####----------------------------------------------------------------------------------------------------------------------------------

if (calculate_repro_number){
  print("Converting to reproductive number")
  repro_spline <- R2spline()
  
  repro_numbers <- lapply(1:nlayers(stacked_pr), function(idx){
    this_layer <- copy(stacked_pr[[idx]])
    this_repro <- apply_spline_to_raster(repro_spline, this_layer, return_r0 = T)
    return(this_repro)
  })
  repro_numbers <- stack(repro_numbers)
  names(repro_numbers) <-  names(stacked_pr)
  
  # Plot and save
  
  # todo: ask Jen for the actual map color bounds
  breaks <- c(0,1,seq(1.5, 5, length.out=25), seq(5.1, 80, length.out=4))
  pal = c("#e0e0e0", terrain.colors(31)[1:30])
  
  if (plot_results){
    pdf(file.path(out_dir, paste0("repro_number_", region, ".pdf")), width=12, height=10)
    
    for(base_val in 1:length(baseline_raster_fnames)){
      end_idx <- (int_count+1) * base_val
      start_idx <- end_idx - int_count
      print(levelplot(repro_numbers[[start_idx:end_idx]], par.settings=rasterTheme(pal), at=breaks,
                      xlab=NULL, ylab=NULL, scales=list(draw=F)) +
              latticeExtra::layer(sp.polygons(africa_shp))
      )
    }
    graphics.off()
    
  }
  
  
  writeRaster(repro_numbers, options="INTERLEAVE=BAND", bylayer=F, suffix="names", filename=file.path(out_dir, paste0("repro_number_",region, suffix, ".tif")), overwrite=T)
}




