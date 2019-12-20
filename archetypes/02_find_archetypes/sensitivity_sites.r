## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# sensitivity_sites.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# December 2019
# 
# For a given covariate and continent of interest, this script rotates the covariate values into n-dimensional space
# as defined by the SVD analysis run in "svd.r" (n is user-defined), then runs k-means in this
# lower-dimensional space and plots the results. 
# 
# For a detailed project write-up see
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(rasterVis)
library(data.table)
library(stats)
library(gridExtra)
library(latticeExtra)
library(Hmisc)
library(FNN)

rm(list=ls())

n_sensitivity_sites <- 10
out_subdir <- "v4_era5_bounded_transmission"

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/")

in_dir <- file.path(base_dir, "results", out_subdir)
guide <- fread(file.path(in_dir, "instructions.csv"))

cluster_df <- unique(guide[, list(continent, n_clusters)])

for (idx in 1:nrow(cluster_df)){
  this_df <- cluster_df[idx]
  
  nclust <- this_df$n_clusters
  # load in relevant k-means outputs
  rotation <- fread(file.path(in_dir, this_df$continent, "02_kmeans", "svd_rotations.csv"))
  load(file.path(in_dir, this_df$continent, "02_kmeans", paste0("k_out_", nclust, "_cluster.rdata")))
  rotation[, cluster:= k_out$cluster]
  cluster_raster <- raster(file.path(in_dir, this_df$continent, "02_kmeans", paste0("map_", nclust, "_cluster", ".tif")))
  
  # find random locations within each cluster
  site_ids <- lapply(1:nrow(k_out$centers), function(this_clust){
    set.seed(206)
    subset <- rotation[cluster==this_clust]
    new_points <- sample(1:nrow(subset), n_sensitivity_sites)
    subset <- subset[new_points]
    return(subset)
  })
  site_ids <- rbindlist(site_ids)
  
  site_id_lat_longs <- data.table(xyFromCell(cluster_raster, site_ids$id))
  setnames(site_id_lat_longs, c("x", "y"), c("longitude", "latitude"))
  site_ids <- cbind(site_ids, site_id_lat_longs)
  
  write.csv(site_ids, file.path(in_dir, this_df$continent, "02_kmeans",  paste0("sensitivity_site_ids_", nclust, "_cluster", ".csv")), row.names=F)
  
  # plot sites on map
  site_id_spoints <- xyFromCell(cluster_raster, site_ids$id, spatial=T)
  
  print("making map")
  palette <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B33539", "#367A40")
  if (nclust<=length(palette)){
    these_colors <- palette[1:nclust]
  }else{
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    these_colors <- gg_color_hue(nclust)
  }
  cluster_raster <- ratify(cluster_raster)
  map_plot <- levelplot(cluster_raster, att="ID", col.regions=these_colors,
                        xlab=NULL, ylab=NULL, scales=list(draw=F),
                        main = "", colorkey=F, margin=F)   +
    latticeExtra::layer(sp.points(site_id_spoints), theme = simpleTheme(col = "black",
                                                                        cex=2))
  
  print(map_plot)
}


