###############################################################################################################
## explore_am_itn_tradeoff.r
## Amelia Bertozzi-Villa
## April 2020
## 
## See how high other interventions need to be to make up for lowered antimalarial coverage
##############################################################################################################

library(data.table)
library(ggplot2)
library(pracma)
library(cgam)
library(raster)
library(rasterVis)

rm(list=ls())

### Setup  ----------------------------------------------------------------------------------------------------------------------------------------------------------

print("Setting up")
func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results")
setwd(func_dir)
source("pr_to_r0.r")
source("map_ii_functions.r")

analysis_subdir <- "20200506_reextract_20191009_mega_era5_new_arch"
analysis_metric <- "prev"

base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
archetype_raster_dir <- file.path(base_dir, "../archetypes/results")
main_dir <- file.path(base_dir, analysis_subdir)

suffix <- ""
out_dir <-  file.path(main_dir,"results", "rasters", "base_2020", analysis_metric)
dir.create(out_dir, recursive = T, showWarnings = F)

region <- "Africa"
africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"
main_raster_input_dir <- "/Volumes/GoogleDrive/My Drive/covid_mitigation_project"
pfpr_raster_fname <- file.path(main_raster_input_dir, "pfpr_2020_realizations", "MEAN.PFPR.2020.tif")

act_raster_fname <- file.path(main_raster_input_dir, "effective_treatment_with_an_antimalarial3/2020.effective.treatment.tif")
itn_raster_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20200418_BMGF_ITN_C1.00_R1.00_V2/04_predictions/ITN_2020_use.tif"

# by what percent do you want to *reduce* act coverage?
act_percent_reductions <- c(20, 40, 60, 80)

### Shapefile  -----------------------------------------------------
print("Loading Shapefile")
africa_shp <- readOGR(africa_shp_dir)
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)

### Desired Geography  -----------------------------------------------------

site_details <- fread(file.path(main_dir, "input", "site_details.csv"))
if (region!="Global"){
  site_details <- site_details[continent==region]
}

# read in cluster rasters, 
archetype_raster <- raster(file.path(archetype_raster_dir, unique(site_details$source), 
                                tolower(region), "02_kmeans",
                                paste0("map_", nrow(site_details),"_cluster.tif")))


### Load and format simulation results  #####----------------------------------------------------------------------------------------------------------------------------------
print("Loading and formatting simulation results")

sim_types <- "prev"
if (analysis_metric %like% "inc"){
  sim_types <- c(sim_types, "inc")
}

sim_result_list <- lapply(sim_types, function(sim_type){
  print(sim_type)
  sim_results <- fread(file.path(main_dir, "results/clean/", paste0("summary_impact_", sim_type, ".csv")))
  sim_results <- sim_results[Site_Name<11]
  
  # subset to just interventions of interest: no new interventions
  sim_results <- sim_results[int_id<=125]
  
  # to start: no irs
  sim_results <- sim_results[irs==0]
  sim_results[, transmission_intensity:=round(log10(x_Temporary_Larval_Habitat), 2)]
  
  sim_results <- unique(sim_results[, list(Site_Name, transmission_intensity, int_id, metric, itn, irs, al_cm, value=mean_final)])
  # merge on a "cluster id" corresponding to the cluster values in each raster.
  # For procedurally generated sites, this should be identical to Site_Name
  # this will also drop the sites that are outside of the region of interest.
  sim_results <- merge(sim_results, site_details[, list(Site_Name=id, cluster_id=cluster)], by="Site_Name")
  return(sim_results)
})
names(sim_result_list) <- sim_types

unique_itn <- sort(unique(sim_result_list[[1]]$itn))
unique_act <- sort(unique(sim_result_list[[1]]$al_cm))

### Load and format rasters  #####----------------------------------------------------------------------------------------------------------------------------------
print("Loading and formatting rasters")
pfpr_raster <- raster(pfpr_raster_fname)
pfpr_raster[is.na(pfpr_raster)] <- 0
itn_raster <- raster(itn_raster_fname)
act_raster <- raster(act_raster_fname)

# crop to itn layer, it's the most restricted
itn_raster[itn_raster==0] <- NA
act_raster <- raster::mask(act_raster, itn_raster)
act_raster[is.na(act_raster) & itn_raster>0] <- 0 #fill in stray NA's

pfpr_raster <- raster::mask(pfpr_raster, itn_raster)

archetype_raster <- extend(crop(archetype_raster, itn_raster), itn_raster)
archetype_raster <- raster::mask(archetype_raster, itn_raster)

print("Converting to data tables")
# convert to data.table for analysis
archetype_dt <- data.table(cellnumber=1:ncell(archetype_raster), archetype=extract(archetype_raster, 1:ncell(archetype_raster)))
archetype_dt <- archetype_dt[!is.na(archetype)]

### Functions  ----------------------------------------------------------------------------------------------------------------------------------------------------------

make_interp_matrix <- function(for_interp, itn_vals, act_vals){
  # for_interp <- this_impact[transmission_intensity==this_intensity]
  # row is itn, column is act
  for_interp <- as.matrix(dcast.data.table(for_interp, itn~al_cm, value.var = "value")[, 2:(length(act_vals)+1)])
  rownames(for_interp) <- itn_vals
  return(for_interp)
}

# for testing
# orig_itn_cov <- 0.4809626
# orig_act_cov <- 0.448775
# new_act <- 0.35902
# prevalence <- 0.7248536
# archetype <- 4
# interp_matrices <- interp_prev_matrices
# interp_matrices_new <- interp_inc_matrices
# act_vals <- unique_act
# itn_vals <- unique_itn


find_new_itn <- function(orig_itn_cov, orig_act_cov, new_act, prevalence, archetype, interp_matrices, interp_matrices_new=NULL,
                         new_act_type="cov", itn_vals=unique_itn, act_vals=unique_act){
  
  if (orig_itn_cov>max(itn_vals)){ # when true coverage is already higher than the max we modeled
    return(-3)
  }

  estimated_prevs <- sapply(interp_matrices[[archetype]], function(this_matrix){
    return(interp2(act_vals, itn_vals, this_matrix, orig_act_cov, orig_itn_cov))
  }) 
  prev_diffs <- abs(estimated_prevs-prevalence)
  intensity_id <- names(which(prev_diffs==min(prev_diffs)))
  # if multiple, pick the highest intensity available
  intensity_id <- intensity_id[length(intensity_id)]
  
  if (intensity_id=="2.2" & prevalence>max(estimated_prevs)){
    return(-4) # simulations cannot replicate this prevalence at this coverage level
  }
  
  if (is.null(interp_matrices_new)){
    val_to_match <- prevalence
    interp_matrices_new <- interp_matrices
  }else{
    val_to_match <- interp2(act_vals, itn_vals, interp_matrices_new[[archetype]][[intensity_id]], orig_act_cov, orig_itn_cov)
  }
  
  new_act <- ifelse(new_act_type=="perc", orig_act_cov*new_act, new_act)
  possible_new_vals <- interp2(act_vals, itn_vals, interp_matrices_new[[archetype]][[intensity_id]], 
                            rep(new_act, length(itn_vals)), itn_vals)
  names(possible_new_vals) <- itn_vals
  possible_new_vals[possible_new_vals<1e-3] <- 0
  
  if (length(unique(possible_new_vals))==1 & unique(possible_new_vals)[[1]]==0){
    this_new_itn_coverage <- -1
  }else if (val_to_match<min(possible_new_vals)){ # when you can't attain the desired prevalence even with 80% itn cov
    this_new_itn_coverage <- -2
  }else if (val_to_match>max(possible_new_vals)){ # when no ITN coverage will still get you less than the target prevalence
    this_new_itn_coverage <- 0
  }else{
    this_new_itn_coverage <- interp1(rev(possible_new_vals), rev(itn_vals), val_to_match)
  }
  return(this_new_itn_coverage)
}

### Main loop: find new ITN coverages  #####----------------------------------------------------------------------------------------------------------------------------------
print("Finding interpolated prevalence matrices")
sim_prevs <- sim_result_list[["prev"]]
interp_prev_matrices <- lapply(unique(sim_prevs$Site_Name), function(this_site){
  print(paste("Site:", this_site))
  these_matrices <- lapply(unique(sim_prevs$transmission_intensity), function(this_intensity){
    for_interp <- sim_prevs[Site_Name==this_site & transmission_intensity==this_intensity]
    return(make_interp_matrix(for_interp, itn_vals = unique_itn, act_vals = unique_act))
  })
  names(these_matrices) <- unique(sim_prevs$transmission_intensity)
  return(these_matrices)
})
names(interp_prev_matrices) <- unique(sim_prevs$Site_Name)

if (analysis_metric %like% "inc"){
  
  print("Finding interpolated incidence matrices")
  sim_incs <- sim_result_list[["inc"]][metric==analysis_metric]
  
  if (analysis_metric %like% "severe"){
    sim_incs[, value:=value*1000]
  }
  
  print("Ensuring monotonicity")
  sim_incs <- rbindlist(lapply(unique(sim_incs$Site_Name), function(this_site){
    this_intensity_incs <- rbindlist(lapply(unique(sim_incs$transmission_intensity), function(this_intensity){
      for_interp <- rbindlist(lapply(unique_act, function(this_cov){
        this_df <- sim_incs[Site_Name==this_site & transmission_intensity==this_intensity & al_cm==this_cov]
        this_df[, new_value:=predict(cgam(value ~ 1+s.decr(itn, numknots = 3), data=this_df, family=gaussian))]
        return(this_df)
      }))
      return(for_interp)
    }))
    return(this_intensity_incs)
  }))
  
  
  pdf(file.path(out_dir, "inc_splines.pdf"), width=28, height=10)
  for (this_site in unique(sim_incs$Site_Name)){
    spline_plot <- ggplot(sim_incs[Site_Name==this_site], aes(x=itn, color=factor(al_cm))) +
                        geom_point(aes(y=value)) +
                        geom_line(aes(y=new_value)) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle=45, hjust=1))+ 
                        facet_grid(al_cm ~ transmission_intensity) +
                        labs(x="ITN Coverage", 
                             y="u10 Severe Incidence (per 1000)",
                             title=paste("Site:", this_site))
    print(spline_plot)
  }
  graphics.off()
  
  setnames(sim_incs, c("value", "new_value"), c("old_value", "value"))
  
  print("Interpolating")
  interp_inc_matrices <- lapply(unique(sim_incs$Site_Name), function(this_site){
    print(paste("Site:", this_site))
    these_matrices <- lapply(unique(sim_incs$transmission_intensity), function(this_intensity){
      for_interp <- sim_incs[Site_Name==this_site & transmission_intensity==this_intensity]
      return(make_interp_matrix(for_interp, itn_vals = unique_itn, act_vals = unique_act))
    })
    names(these_matrices) <- unique(sim_incs$transmission_intensity)
    return(these_matrices)
  })
  names(interp_inc_matrices) <- unique(sim_incs$Site_Name)
}

### Main loop: find new ITN coverages  #####----------------------------------------------------------------------------------------------------------------------------------

raster_template <- copy(archetype_raster)
raster_template[raster_template>0] <- NA

print("Converting coverage and pfpr to data tables")
raster_stack <- stack(pfpr_raster, itn_raster, act_raster)
names(raster_stack) <- c("prevalence", "itn_orig", "act_orig")
raster_dt <- cbind(archetype_dt, extract(raster_stack, archetype_dt$cellnumber))


# note: these take about an hour each, be warned
for(percent_reduction in act_percent_reductions){
  print(paste("Reducing act coverage by", percent_reduction))
  new_act_prop <- 1-(percent_reduction/100)
  
  this_raster_dt <- copy(raster_dt)
  this_raster_dt[, act_prop:=new_act_prop]
  this_raster_dt[, act_new:=act_orig*new_act_prop]
  print("Calculating new ITN coverage")
  tic <- Sys.time()
  print(tic)
  
  # for (idx in 1:nrow(this_raster_dt)){
  #   print(idx)
  #   this_dt <- this_raster_dt[idx]
  #   this_dt[, itn_new:=find_new_itn(itn_orig, act_orig, act_new, prevalence, archetype, 
  #                                          interp_matrices = interp_prev_matrices,
  #                                          interp_matrices_new=interp_inc_matrices
  #   ),
  #   by = seq_len(nrow(this_dt))]
  # }
  # 
  
  
  if (analysis_metric %like% "inc"){
    this_raster_dt[, itn_new:=find_new_itn(itn_orig, act_orig, act_new, prevalence, archetype, 
                                           interp_matrices = interp_prev_matrices,
                                           interp_matrices_new=interp_inc_matrices
                                           ),
                   by = seq_len(nrow(this_raster_dt))]
  }else{
    this_raster_dt[, itn_new:=find_new_itn(itn_orig, act_orig, act_new, prevalence, archetype, interp_matrices = interp_prev_matrices),
                   by = seq_len(nrow(this_raster_dt))]
  }
  
  toc <- Sys.time()
  elapsed <- toc-tic
  print(paste("Time elapsed:", elapsed, units(elapsed)))
  
  this_raster_dt[itn_new>0, itn_final:=itn_new]
  this_raster_dt[itn_new==0, itn_final:=itn_orig] # better alternative than setting to zero
  this_raster_dt[itn_new==-1, itn_final:=itn_orig] # places where we project no transmission (archetype 1)
  this_raster_dt[itn_new==-2, itn_final:=1] # maximize nets since we don't know the upper bound
  this_raster_dt[itn_new==-3, itn_final:=itn_orig] # places where cov is already >80%
  this_raster_dt[itn_new==-4, itn_final:=1] # placeholder until we get a better idea

  
  write.csv(this_raster_dt, file.path(out_dir, paste0("derived_itn_cov_reduce_act_", percent_reduction, "_", analysis_metric, ".csv")), row.names = F)
  
  this_itn_raster <- copy(raster_template)
  this_itn_raster[this_raster_dt$cellnumber] <- this_raster_dt$itn_final
  this_itn_raster[is.na(this_itn_raster) & !is.na(itn_raster)] <- itn_raster[is.na(this_itn_raster) & !is.na(itn_raster)]
  
  writeRaster(this_itn_raster, filename = file.path(out_dir, paste0("derived_itn_cov_reduce_act_", percent_reduction, "_", analysis_metric, ".tif")), overwrite=T)
  
}




# # pdf("~/Desktop/itn_am_heatmaps.pdf", width=14, height=8)
# 
# for (this_site in unique(impact$Site_Name)){
#   this_heatmap <- 
    # ggplot(sim_prevs[Site_Name==4 & transmission_intensity>-4], aes(x=al_cm, y=itn, fill=value)) +
    # geom_tile() +
    # geom_text(aes(label=round(value, 2))) +
    # scale_fill_distiller(name="Incidence", palette = "Spectral") +
    # facet_wrap(.~transmission_intensity) +
    # theme_minimal() +
    # scale_x_continuous(breaks=seq(0, 0.8, 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) +
    # scale_y_continuous(breaks=seq(0, 0.8, 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) +
    # labs(x="Effective Treatment",
    #      y="ITN Coverage")
#   # print(this_heatmap)
# }
# 
# # graphics.off()








