## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 02_k_means.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# May 2018
# 
# For a given covariate and continent of interest, this script rotates the covariate values into n-dimensional space
# as defined by the SVD analysis run in "svd.r" (n is user-defined), then runs k-means in this
# lower-dimensional space and plots the results. 
# 
# For a detailed project write-up see
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(data.table)
library(stats)
library(gridExtra)
library(latticeExtra)
library(Hmisc)
library(FNN)

set.seed(206)

rm(list=ls())
overwrite_rotation <- F
overwrite_kmeans <- T
out_subdir <- "original_megatrends"

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/")

out_dir <- file.path(base_dir, "results", out_subdir)
guide <- fread(file.path(out_dir, "instructions.csv"))

cluster_counts <- 3:20

for (this_continent in guide$continent){
  
  this_cov <- gsub("/", ".", guide[continent==this_continent]$covariates)
  print(paste("clustering", this_cov, "for", this_continent))
  
  nvecs <- guide[continent==this_continent]$singular_vectors
  
  this_out_dir <- file.path(out_dir, this_continent, "02_kmeans")
  dir.create(this_out_dir, showWarnings=F, recursive=T)
  svd_dir <- file.path(this_out_dir, "../01_svd")
  cov_dir <- file.path(base_dir, "covariates", unique(guide$cov_directory), this_continent)
  
  # rotate matrix, if needed
  rotation_fname <- file.path(this_out_dir, "svd_rotations.csv")
  
  if (file.exists(rotation_fname) & overwrite_rotation==F){
    print("loading matrix rotations")
    rotation <- fread(rotation_fname)
    # loads svd outputs to get dataset of cell values
    load(file.path(svd_dir, "svd_output.rdata")) 
  }else{
    print("finding matrix rotations")
    
    # loads svd_out (svd outputs, including singular vectors) and svd_wide_datatable (the data used to run svd)
    load(file.path(svd_dir, "svd_output.rdata")) 
    sing_vecs <- svd_out$u[, 1:nvecs]
    
    ## multiply by original matrix to get rotations
    print("rotating")
    rotation <- data.frame(t(t(sing_vecs)%*%as.matrix(svd_wide_datatable[,4:ncol(svd_wide_datatable)])))
    rotation$id <- as.integer(rownames(rotation))
    rotation <- data.table(rotation)
    write.csv(rotation, rotation_fname, row.names=F)
    
  }
  
  all_vals <- melt(svd_wide_datatable, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
  all_vals[, id:=as.integer(as.character(id))]
  rm(svd_out, svd_wide_datatable); gc()
  
  # find map from raster cells to rotation ids
  temp_raster <- raster(file.path(cov_dir, "mask.tif"))
  
  # bit of matrix transposition to appropriately map a raster ID to a matrix cell
  print("finding raster cell: rotation id map")
  for_centroids <- rep(NA, ncell(temp_raster))
  for_centroids[rotation$id] <- rotation$id
  for_centroids <- matrix(for_centroids, nrow=nrow(temp_raster))
  for_centroids <- base::t(for_centroids)
  cell_id_map <- data.table(raster_cell = which(!is.na(for_centroids)),
                            id = for_centroids[!is.na(for_centroids)])
  write.csv(cell_id_map, file.path(this_out_dir, "cell_id_map.csv"), row.names=F)
  
  
  # k-means
  for (nclust in cluster_counts){
    
    # if k-means has already been run, just load outputs
    k_out_fname <- file.path(this_out_dir, paste0("k_out_", nclust, "_cluster", ".rdata"))
    cluster_raster_fname <- file.path(this_out_dir, paste0("map_", nclust, "_cluster", ".tif"))
    random_trace_fname <- file.path(this_out_dir, paste0("random_trace_", nclust, "_cluster",  ".csv"))
    summary_fname <- file.path(this_out_dir,  paste0("summary_", nclust, "_cluster", ".csv"))
    site_id_fname <- file.path(this_out_dir,  paste0("site_ids_", nclust, "_cluster", ".csv"))
    
    if (file.exists(k_out_fname) & overwrite_kmeans==F){
      print(paste("k-means already run for cluster count", nclust))
    }else{
      print(paste("finding clusters for k of", nclust))
      k_out <- kmeans(rotation[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
      rotation[, cluster:= k_out$cluster]
      
      print("creating new raster")
      cluster_raster <- rep(NA, ncell(temp_raster))
      cluster_raster[rotation$id] <- rotation$cluster
      cluster_raster <- matrix(cluster_raster, nrow=nrow(temp_raster))
      cluster_raster <- raster(cluster_raster, template=temp_raster)
      writeRaster(cluster_raster, cluster_raster_fname, overwrite=T)
  
      all_vals <- merge(all_vals, rotation[, list(id, cluster)], by="id", all=T)
      
      print("finding random traces")
      random_grids <- sample(unique(all_vals$id), 500)
      random_trace <- all_vals[id %in% random_grids]
      write.csv(random_trace, random_trace_fname, row.names=F)
    
      print("finding summary stats")
      summary_vals <- all_vals[,list(mean=mean(cov_val),
                                     perc_25=quantile(cov_val, c(0.25)),
                                     perc_75=quantile(cov_val, c(0.75)),
                                     perc_05=quantile(cov_val, c(0.025)),
                                     perc_95=quantile(cov_val, c(0.975)),
                                     median=median(cov_val)), by=list(cluster, variable_name, variable_val, cov)]
      summary_vals[, nclust:=nclust]
      write.csv(summary_vals, summary_fname, row.names=F)
      
      # delete so all_vals can be used again for the next round of k-means
      all_vals[, cluster:=NULL]
      
      print("finding centroids for representative sites")
      # find cell IDs of centroids and map them to the appropriate locations
      site_ids <- rotation[get.knnx(rotation[,1:nvecs], k_out$centers, k=1)$nn.index]
      site_ids <- merge(site_ids, cell_id_map, by="id", all.x=T)
      write.csv(site_ids, site_id_fname, row.names=F)
      
      save(k_out, file=k_out_fname) # save last to ensure other outputs have been saved if this exists
      
    }
    
  }

}


