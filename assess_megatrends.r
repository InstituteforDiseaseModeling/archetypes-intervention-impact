## -----------------------------------------------------------------------------------------------------------------
# Megatrends: Characterizing residual transmission
# assess_megatrends.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# December 2018
# 
# As per the request of the WHO-SAGme, look at the megatrends + existing intervention forecasts and
# take a first pass at determining:
#   1. Whether pockets of residual transmission are determined by high initial transmission or something else;
#   2. If "something else", what mechanistic modeling scenarios result in higher residual transmission.
# 
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(rasterVis)

rm(list=ls())

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact")

main_dir <- file.path(base_dir, "writing_and_presentations/megatrends/pfpr_rasters")
cluster_fname <- file.path(base_dir, "lookup_tables/interactions/africa_clusters_v4.tif")


# load in data, clip to africa
cluster_layer <- raster(cluster_fname)
megatrends_only <- raster(file.path(main_dir, "actual_ssp2_noint_2050.tif"))
megatrends_full <- raster(file.path(main_dir, "actual_ssp2_2050.tif"))
# interventions <- raster(file.path(main_dir, "actual_ssp2_2050_ITN80ACT80-14.tif"))
interventions <- raster(file.path(main_dir, "pfpr_africa.tif"), band=2)
cluster_layer <- crop(cluster_layer, interventions)
megatrends_only <- crop(megatrends_only, cluster_layer)
megatrends_full <- crop(megatrends_full, cluster_layer)
interventions <- crop(interventions, cluster_layer)

# mask clusters to megatrends
cluster_layer <- raster::mask(cluster_layer, megatrends_only)

cutoff_pr <- 0.001

# mask areas that are zero in megatrends only, these will by definition be zero in intervenions as well
megatrends_orig <- copy(megatrends_only)
megatrends_only[megatrends_orig<cutoff_pr] <- -Inf
megatrends_full[megatrends_orig<cutoff_pr] <- -Inf
interventions[megatrends_orig<cutoff_pr] <- -Inf

bounded_interventions <-min(stack(interventions,megatrends_full))
megatrends_diff <- interventions - bounded_interventions

# explore areas of residual transmission--------------------------

resid_interventions <- copy(interventions)
resid_megatrends <- copy(megatrends_only)
resid_clusters <- copy(cluster_layer)

resid_interventions[interventions<cutoff_pr] <- -Inf
resid_megatrends[interventions<cutoff_pr] <- -Inf
resid_clusters[interventions<cutoff_pr] <- -Inf
writeRaster(resid_clusters, file.path(main_dir, "masked_clusters.tif"), overwrite=T)

reduction <- resid_megatrends-resid_interventions

raster_to_dt <- function(vals, name){
  vals <- as.matrix(vals)
  vals <- data.table(id = which(!is.na(vals) & !is.infinite(vals)),
                     value = vals[!is.na(vals) & !is.infinite(vals)])
  setnames(vals, "value", paste0(name, "_val"))
  return(vals)
}

resid_megatrends_dt <- raster_to_dt(resid_megatrends, name="megatrend")
resid_interventions_dt <- raster_to_dt(resid_interventions, name="intervention")
resid_cluster_dt <- raster_to_dt(resid_clusters, name="cluster")

reduction_dt <- merge(resid_megatrends_dt, resid_interventions_dt, by="id", all=T)
reduction_dt <- merge(reduction_dt, resid_cluster_dt, by="id")
reduction_dt[, cluster_val:=as.factor(cluster_val)]
reduction_dt[, perc_reduction:= (megatrend_val-intervention_val)/megatrend_val *100]

colors <- c("#00a08a", "#d71b5a", "#f2a200", "#f98400", "#902e57", "#5392c2")

ggplot(reduction_dt, aes(x=megatrend_val, y=intervention_val)) + 
        geom_point(aes(color=cluster_val), alpha=0.5) + 
        scale_color_manual(values=colors) + 
         facet_wrap(~cluster_val) + 
        geom_abline() +
        theme_minimal() + 
        theme(legend.position="none") + 
        labs(title="Megatrend Only vs Megatrend Plus Intervention PfPR, \n by Transmission Archetype",
             x="Megatrend Only",
             y="Megatrend + ITN 80%, ACT 80%")


ggplot(reduction_dt, aes(x=megatrend_val, y=perc_reduction)) + 
  geom_point(aes(color=cluster_val), alpha=0.5) + 
  scale_color_manual(values=colors) + 
  facet_wrap(~cluster_val) + 
  theme_minimal() + 
  theme(legend.position="none") + 
  labs(title="Megatrend Only PfPR vs % Reduction under ITN 80%, ACT 80%, \n by Transmission Archetype",
       x="Megatrend Only",
       y="PFPr % Reduction under ITN 80%, ACT 80%")

# todo: spatially disaggregate pixels that track lookup table vs megatrends 

## todo: better map of where residual transmission is
resid_clusters <- raster(file.path(main_dir, "masked_clusters.tif"))
resid_clusters <- ratify(resid_clusters)
levelplot(resid_clusters, att="ID", col.regions=colors,
          xlab=NULL, ylab=NULL, scales=list(draw=F),
          main = "", colorkey=F, margin=F)
