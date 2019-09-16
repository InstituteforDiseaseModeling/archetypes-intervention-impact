library(gdistance)
library(data.table)
library(stringr)
library(stats)
library(rasterVis)
library(ggplot2)
library(gridExtra)

rm(list=ls())

palette <- "Paired"
source("classify_functions.r")
base_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification",
                      "vectors")
mask_dir <- "Z:/mastergrids/Global_Masks/MAP_Regions/MAP_Regions_Pf_5k.tif"
vector_dir <- "Z:/master_publication_datasets/SinkaEtAl_AfricanVectorRelativeAbundance_MalJ_2016/"

in_fnames <- c("arabiensis", "gambiae", "funestus")
mask_raster <- get_mask("africa", file.path(base_dir), mask_dir)

vector_list <- NULL

# extract rasters
for (idx in 1:3){
  
  species <- in_fnames[idx]
  print(paste("clipping global", species, "raster"))
  out_fname <- file.path(base_dir, paste0(species, ".tif"))
  
  if (!file.exists(out_fname)){
    full <- raster(file.path(vector_dir, "All-multiband.RelativeAbundance.tif"), band=idx)
    vals <- crop(full, mask_raster)
    mask_raster <- crop(mask_raster, vals) # ensure that extents of two rasters are the same
    vals <- mask(vals, mask_raster, maskvalue=FALSE)
    writeRaster(vals, out_fname, overwrite=T)
  }else{
    vals <- raster(out_fname)
  }
  
  vals <- as.matrix(vals)
  vals <- data.table(species=species, 
                     id = which(!is.na(vals)),
                     abundance = vals[!is.na(vals)])
  
  vector_list[[idx]] <- vals
}

all_vectors <- rbindlist(vector_list)

# remove pixels with no mosquitoes
all_vectors[, tot:=sum(abundance), by="id"]
all_vectors <- all_vectors[tot!=0]
all_vectors[, tot:=NULL]

# k-means 
nclust <- 4

pdf(file.path(base_dir, "k_means_vectors.pdf"), width=9, height=6)
for (nclust in 3:7){
  print(paste("finding clusters for k of", nclust))
  for_kmeans <- dcast(all_vectors, id~species, value.var="abundance")
  k_out <- kmeans(for_kmeans[, 2:4], centers=nclust, algorithm = "MacQueen", iter.max=100)
  for_kmeans[, cluster:= k_out$cluster]
  
  print("creating new raster")
  # load mask raster to get dimensions & extent
  cluster_raster <- matrix(nrow=mask_raster@nrows, ncol=mask_raster@ncols)
  cluster_raster[for_kmeans$id] <- for_kmeans$cluster
  cluster_raster <- raster(cluster_raster, template=mask_raster)
  
  cluster_raster <- ratify(cluster_raster)
  map_plot <- levelplot(cluster_raster, att="ID", col.regions=brewer.pal(nclust, palette),
                        xlab=NULL, ylab=NULL, scales=list(draw=F),
                        main = paste("Mapped Clusters"),  margin=F)
  
  for_bar <- melt(for_kmeans, id.vars = c("id", "cluster"), variable.name="species", value.name="rel_abundance")
  for_bar <- for_bar[, list(rel_abundance=mean(rel_abundance)), by=c("cluster", "species")]
  
  bar_plot <- ggplot(for_bar, aes(x=cluster, y=rel_abundance)) +
                geom_bar(aes(fill=species), stat="identity") + 
                theme_minimal() +
                theme(legend.position="bottom") +
                scale_fill_manual(values=c("#ffd92f", "#a6d854", "#8da0cb"), name="")+ 
                labs(title="Relative Species Abundance",
                     x="Cluster",
                     y="")
  
  
  layout <- rbind(c(1,1,NA),
                  c(1,1,2),
                  c(1,1,2),
                  c(1,1,NA))
  
  print("saving")
  full_plot <- grid.arrange(map_plot, bar_plot, layout_matrix=layout)
  print(full_plot)
}
graphics.off()

# writeRaster(cluster_raster, cluster_raster_fname, overwrite=T)








