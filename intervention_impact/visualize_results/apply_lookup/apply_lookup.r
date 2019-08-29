
library(data.table)
library(raster)
library(rasterVis)

rm(list=ls())

# todo: more explicit file paths, take things out of this repo folder where they don't belong

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results/apply_lookup")
setwd(func_dir)
source("pr_to_r0.r")
source("apply_lookup_functions.r")
out_dir <- "gates_output"

### Input variables  -----------------------------------------------------

calculate_par <- F
pop_fname <- "other/ssp5_total_2050_MG_5K.tif"

calculate_repro_number <- T 

region <- "africa" # geography for which you want to predict. currently responds only to "africa" or "global"

# what interventions from the lookup table do you want to run?
interventions <- c("ITN 0.8; ACT 0.8;",
                   "ITN 0.8; IRS 0.8; ACT 0.8;") 

# other intvervention values I need:
# Maximized current tools + ATSB 0.15%, 3%: lookup_outdoor_interventions_sweep_timing
# Maximized current tools + ATSB 15%, 25%:  atsb_to_eliminate.csv
# Drugs & vaccines: lookup_full_interactions_experimental.csv
# ATSB 1% at different baselines: see other/apply_lookup_scratch.r

# what rasters do you want to use as baselines? You'll get a new set of results for each
baseline_raster_fnames <- list( # "True PfPR 2017"="PfPR_rmean_Global_admin0_2017.tif",
  "Megatrends Base 2000"="actual_ssp2_base2000_2050.tif",
  "Megatrends Base 2016"="actual_ssp2_base2016_2050.tif"
)

# What raster do you want to use as a maximum value on all results?
bounding_fname <- "actual_ssp2_base2016_2050.tif"

# Are there additional rasters you want to visualize, but not apply a lookup table to?
comparison_fnames <- c()

### Calculate reproductive number, or just PR?  -----------------------------------------------------


### Desired Geography  -----------------------------------------------------

# read in cluster rasters, assign values from LUT
cluster_map <- raster(paste0(region, "_clusters_v4.tif"))

# map of cluster numbers to site names in the lut
cluster_list <-  c("aba"=1,
                      "kasama"=3,
                      "djibo"=4,
                      "kananga"=2,
                      "moine"=6,
                      "gode"=5)


if (region=="global"){
  cluster_list <- c(cluster_list, c("bajonapo"=7, "karen"=8))
}

cluster_key <- data.table(Site_Name=names(cluster_list),
                          cluster_id=cluster_list)

### Lookup table -----------------------------------------------------

# read in lookup table
lut <- fread("lookup_full_interactions_v4.csv")

# keep only the desired interventions; simplify
lut <- lut[Intervention %in% interventions]
lut <- unique(lut[, list(Site_Name, Intervention, mean_initial, mean_final)])
lut <- merge(lut, cluster_key, by="Site_Name", all.x=T)

### Bounding raster -----------------------------------------------------
bounding_pr <- raster(bounding_fname)
cluster_map <- crop(cluster_map, bounding_pr)
bounding_pr <- crop(bounding_pr, cluster_map)
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
  baseline_pr <- raster(baseline_raster_fnames[[baseline_label]])
  baseline_pr <- crop(baseline_pr, cluster_map)
  baseline_pr <- raster::mask(baseline_pr, cluster_map)
  names(baseline_pr) <- baseline_label
  
  pr_list[[pr_idx]] <- baseline_pr
  pr_idx <- pr_idx + 1 
  
  for (intervention in interventions){
    print(intervention)
    pr_list[[pr_idx]] <- apply_lookup(intervention, baseline_pr, lut, cluster_map, bounding_pr)
    pr_idx <- pr_idx + 1 
  }
}

stacked_pr <- stack(pr_list)

# Plot and save
color_vals <- generate_full_pal()
row_count <- length(interventions)

# pdf(file.path(out_dir, paste0("pfpr_", region, ".pdf")), width=12, height=6)
print(levelplot(stacked_pr, par.settings=rasterTheme(color_vals$pal), at=color_vals$breaks,
                xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), layout=c(row_count+1, row_count))
)
# graphics.off()
# writeRaster(stacked_pr, options="INTERLEAVE=BAND", bylayer=T, suffix="names", filename=file.path(out_dir, paste0("pfpr_",region, ".tif")), overwrite=T)


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
  pal <- wpal("seaside")
  # pdf(file.path(out_dir, paste0("atsb_par_", continent, ".pdf")), width=12, height=6)
  print(levelplot(par + 1, par.settings=rasterTheme(pal), zscaleLog=T, margin=F, 
                  xlab=NULL, ylab=NULL, scales=list(draw=F), layout=c(row_count+1, row_count))) # +
          # latticeExtra::layer(sp.polygons(africa_shp)))

  # graphics.off()
  
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
  # pdf(file.path(out_dir, paste0("RC_", region, ".pdf")), width=12, height=6)
  print(levelplot(repro_numbers, par.settings=rasterTheme(pal), at=breaks,
                  xlab=NULL, ylab=NULL, scales=list(draw=F), layout=c(row_count+1, row_count))
        )
  # graphics.off()
  # writeRaster(stacked_pr, options="INTERLEAVE=BAND", bylayer=T, suffix="names", filename=file.path(out_dir, paste0("RC_",region, ".tif")), overwrite=T)
}




