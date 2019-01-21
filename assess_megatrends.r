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
library(gridExtra)

rm(list=ls())

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact")

main_dir <- file.path(base_dir, "writing_and_presentations/megatrends/pfpr_rasters")
cluster_fname <- file.path(base_dir, "lookup_tables/interactions/africa_clusters_v4.tif")

colors <- c("#00a08a", "#d71b5a", "#f2a200", "#f98400", "#902e57", "#5392c2")

# load in data, clip to africa
cluster_layer <- raster(cluster_fname)
megatrends_noints <- raster(file.path(main_dir, "actual_ssp2_noint_2050.tif"))
megatrends_ints <- raster(file.path(main_dir, "actual_ssp2_2050.tif"))

int_type <- "itn_act"

if (int_type=="itn_act"){
  interventions <- raster(file.path(main_dir, "pfpr_africa.tif"), band=2)
  interventions_pete <- raster(file.path(main_dir, "actual_ssp2_2050_ITN80ACT80-14.tif"))
  int_label <- "ITN 60%, ACT 60%"
}else if (int_type=="itn_irs_act"){
  interventions <- raster(file.path(main_dir, "pfpr_africa.tif"), band=3)
  interventions_pete <- raster(file.path(main_dir, "actual_ssp2_2050_ITN80IRS80ACT80-14.tif"))
  int_label <- "ITN 60%, IRS 60%, ACT 60%"
}else{
  stop(paste("unrecognized intervention type", int_type))
}

# ensure consistent extents in all rasters, and
# mask areas that are zero (or near-zero) in megatrends. 
# These will by definition be zero in interventions as well

cutoff_pr <- 0.0001
cluster_layer <- crop(cluster_layer, interventions)
megatrends_orig <- copy(megatrends_noints)

rastlist <- list("megatrends_noints"=megatrends_noints, "megatrends_ints"=megatrends_ints,
                 "interventions"=interventions, "interventions_pete"=interventions_pete)

print("cropping and removing zeros")
for (rastname in names(rastlist)){
  rast <- rastlist[[rastname]]
  rast[rast<cutoff_pr] <- 0
  newrast <- crop(rast, cluster_layer)
  # save and reload such that -Inf gets applied appropriately
  # todo: find a better way to reset the "names" variable in raster layers
  fname <- file.path(main_dir, "temp", paste0(rastname, ".tif"))
  writeRaster(newrast, fname, overwrite=T)
  newrast <- raster(fname)
  assign(rastname, newrast)
}

# mask clusters to megatrends
cluster_layer <- raster::mask(cluster_layer, megatrends_noints)
names(cluster_layer) <- "cluster"

# my "interventions" dataset is unbounded by the "megatrends with interventions" dataset.
# I need to bound it to be comparable to Pete's 
# (and to explore the difference between bounded and unbounded residual transmission)

# TODO: find reason for differences between my and Pete's bounded rasters. 
bounded_interventions <-min(stack(interventions,megatrends_ints))
names(bounded_interventions) <- "bounded_interventions"
megatrends_diff <- interventions - bounded_interventions

# compare my bounded intervention to Pete's
abv_pg_diff <- bounded_interventions - interventions_pete
plot(abv_pg_diff>0.01 | abv_pg_diff<(-0.01), main=">1% difference between ABV and PG")


# spatially disaggregate pixels that track lookup table vs megatrends 
pdf(file.path(main_dir, "cluster_plots.pdf"), width=12, height=6)

for (cluster_idx in 1:6){
  this_color <- colors[[cluster_idx]]
  stacked_series <- stack(cluster_layer==cluster_idx,
                          cluster_layer==cluster_idx & bounded_interventions>cutoff_pr,
                          cluster_layer==cluster_idx & bounded_interventions>cutoff_pr & bounded_interventions<interventions)
  
  names(stacked_series) <- c("Full Cluster", "Residual Transmission", "Megatrends RT")
  cluster_plot <- levelplot(stacked_series, att="ID", col.regions=c("#A9A9A9", this_color),
                            xlab=NULL, ylab=NULL, scales=list(draw=F),
                            main = int_label, colorkey=F, margin=F)
  print(cluster_plot)
}
graphics.off()



# explore areas of residual transmission--------------------------


raster_to_dt <- function(rast){
  vals <- as.matrix(rast)
  vals <- data.table(id = which(!is.na(vals) & !is.infinite(vals)),
                     value = vals[!is.na(vals) & !is.infinite(vals)])
  vals[, type:= names(rast)]
  return(vals)
}

reduction_dt <- lapply(list(interventions, bounded_interventions, megatrends_ints, megatrends_noints, cluster_layer), raster_to_dt)
reduction_dt <- rbindlist(reduction_dt)
reduction_dt <- dcast.data.table(reduction_dt, id ~ type)
reduction_dt <- reduction_dt[complete.cases(reduction_dt)] # todo: what are the 10,000 pixels that get dropped here between bounded_interventions and interventions?

reduction_dt[, cluster:=as.factor(cluster)]
reduction_dt[, perc_reduction:= (megatrends_noints-bounded_interventions)/megatrends_noints *100]

ggplot(reduction_dt, aes(x=megatrends_noints, y=megatrends_ints)) +
  geom_point(alpha=0.25) +
  facet_wrap(~cluster) +
  geom_abline() +
  theme_minimal() +
  theme(legend.position="none") +
  xlim(0,1)  +
  labs(title="Megatrend (Base 2000) vs Megatrend (Base 2016), \n by Transmission Archetype",
       x="Megatrends, Base 2000",
       y="Megatrends, Base 2016")

ggplot(reduction_dt[interventions>cutoff_pr], aes(x=megatrends_noints, y=bounded_interventions)) + 
        geom_point(aes(y=megatrends_ints), alpha=0.1) + 
        geom_point(aes(color=cluster), alpha=0.5) + 
        scale_color_manual(values=colors) + 
        facet_wrap(~cluster) + 
        geom_abline() +
        theme_minimal() + 
        theme(legend.position="none") + 
        xlim(0,1)  + 
        labs(title="Megatrend (Base 2000) vs Megatrend Plus Interventions, \n by Transmission Archetype",
             x="Megatrend, Base 2000",
             y=paste("Megatrend Base 2000 +", int_label))


## some putzing around --------------------------
reduction_dt[, type:=ifelse(interventions==bounded_interventions, "lookup_val", "megatrends_val")]

summary_pulls <- reduction_dt[, .N, by=list(cluster, type)]
summary_pulls <- dcast(summary_pulls, cluster~type)
summary_pulls[, perc_lookup:=(lookup_val/(megatrends_val+lookup_val))*100]

test <- reduction_dt[cluster==5 & bounded_interventions > 0.001]
test[, id:=factor(id)]
ggplot(test, aes(x=id)) +
  geom_point(aes(y=megatrends_noints), color="red") +
  geom_point(aes(y=bounded_interventions), color="blue") 
  # geom_point(aes(y=interventions), color="green") 

## ------------------------------------------------------


