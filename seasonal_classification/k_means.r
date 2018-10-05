## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# k_means.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# May 2018
# 
# For a given covariate and continent of interest, this script rotates the covariate values into n-dimensional space
# as defined by the SVD analysis run in "extract-and_svd.r" (n is user-defined), then runs k-means in this
# lower-dimensional space and plots the results. 
# 
# For a detailed project write-up see
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------


library(rasterVis)
library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(gridExtra)
library(Hmisc)

set.seed(206)

rm(list=ls())

source("classify_functions.r")
base_dir <- file.path(Sys.getenv("USERPROFILE"), "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")
palette <- "Paired" # color scheme for plots 

# number of singular vectors to use, from visual inspection of svd plots
# input_list <- list(tsi= list(africa=2, asia=3, americas=2),
#                    rainfall = list(africa=3, asia=3, americas=3)) 

input_list <- list(tsi_rainfall_vector_abundance=list(africa=3))

for (this_cov in names(input_list)){
  nvec_list <- input_list[[this_cov]]
  
  for (continent in names(nvec_list)){
    print(paste("clustering", this_cov, "for", continent))
    
    main_dir <- file.path(base_dir, continent)
    
    # rotate matrix, if needed
    rotation_fname <- file.path(main_dir, paste0("svd_rotations_", this_cov, ".csv"))
    nvecs <- nvec_list[[continent]]
    
    if (file.exists(rotation_fname)){
      print("loading matrix rotations")
      rotation <- fread(rotation_fname)
    }else{
      print("finding matrix rotations")
      rotation <- rotate_matrix(nvecs, main_dir, this_cov)
      write.csv(rotation, rotation_fname, row.names=F)
    }
    
    all_vals <- fread(file.path(main_dir, paste0(this_cov, "_vals.csv")))
    setnames(all_vals, "value", "cov_val")
    
    # k-means and plotting
    pdf(file.path(main_dir, paste0("k_means_", this_cov, ".pdf")), width=9, height=6)
    for (nclust in 3:7){
      
      cov_label <- ifelse(nchar(this_cov)==3, toupper(this_cov), capitalize(this_cov))
      
      # if k-means has already been run, just load outputs
      k_out_fname <- file.path(main_dir, paste0("k_out_", this_cov, "_", nclust, ".tif"))
      cluster_raster_fname <- file.path(main_dir, paste0("k_clusters_", this_cov, "_", nclust, ".tif"))
      random_trace_fname <- file.path(main_dir, paste0("random_trace_", this_cov, "_", nclust, ".csv"))
      time_series_fname <- file.path(main_dir,  paste0("time_series_", this_cov, "_", nclust, ".csv"))
      
      if (file.exists(k_out_fname)){
        print("k-means already run, loading outputs")
        cluster_raster <- raster(cluster_raster_fname)
        random_trace <- fread(random_trace_fname)
        time_series <- fread(time_series_fname)
        
      }else{
        print(paste("finding clusters for k of", nclust))
        k_out <- kmeans(rotation[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
        rotation[, cluster:= k_out$cluster]
        
        print("creating new raster")
        # load mask raster to get dimensions & extent
        if (this_cov %like% "vector_abundance"){
          temp_raster <- raster(file.path(main_dir, "rasters", "vector_abundance_species_gambiae.tif"))
        }else{
          temp_raster <- raster(file.path(main_dir, "rasters", "rainfall_month_01.tif"))
        }
        
        cluster_raster <- rep(NA, ncell(temp_raster))
        cluster_raster[rotation$id] <- rotation$cluster
        cluster_raster <- matrix(cluster_raster, nrow=nrow(temp_raster))
        cluster_raster <- raster(cluster_raster, template=temp_raster)
        writeRaster(cluster_raster, cluster_raster_fname, overwrite=T)
    
        all_vals <- merge(all_vals, rotation[, list(id, cluster)], by="id", all=T)
        
        print("finding random traces")
        random_grids <- sample(unique(all_vals$id), 500)
        random_trace <- all_vals[variable_name=="month" & id %in% random_grids]
        random_trace[, variable_val:=as.integer(variable_val)]
        write.csv(random_trace, random_trace_fname, row.names=F)
        
        print("finding mean time series")
        time_series <- all_vals[variable_name=="month", list(cov_val=mean(cov_val)), by=list(cluster, variable_val, cov)]
        time_series[, nclust:=nclust]
        time_series[, variable_val:= as.integer(variable_val)]
        write.csv(time_series, time_series_fname, row.names=F)
        all_vals[, cluster:=NULL]
        
        save(k_out, file=k_out_fname) # save last to ensure other outputs have been saved if this exists
        
      }
      
      plotlist <- NULL
      
      print("making map")
      cluster_raster <- ratify(cluster_raster)
      map_plot <- levelplot(cluster_raster, att="ID", col.regions=brewer.pal(nclust, palette),
                            xlab=NULL, ylab=NULL, scales=list(draw=F),
                            main = paste("Clusters"), colorkey=F, margin=F)
      plotlist[[1]] <- map_plot
      
      lines <- lapply(unique(time_series$cov), function(cov_value){
        ggplot(time_series[cov==cov_value], aes(x=variable_val, y=cov_val)) +
          geom_line(aes(color=factor(cluster)), size=2) +
          geom_line(data=random_trace[cov==cov_value], aes(group=id, color=factor(cluster)), size=0.5, alpha=0.25) + 
          facet_grid(cluster~.) + 
          scale_color_brewer(type="qual", palette=palette) +
          scale_x_continuous(breaks=1:12, labels=1:12) +
          theme_minimal() +
          theme(legend.position="none",
                plot.title = element_text(face="bold")) +
          labs(title=paste(cov_value),
               x="Month",
               y="")
      })
      
      plotlist <- append(plotlist, lines)
      
      layout_pattern <- c(1, 1:(length(lines)+1))
      layout <- rbind(layout_pattern,
                      layout_pattern)
      
      print("saving")
      full_plot <- grid.arrange(grobs=plotlist, layout_matrix=layout)
      print(full_plot)
      
    }
    graphics.off()
  }
}

