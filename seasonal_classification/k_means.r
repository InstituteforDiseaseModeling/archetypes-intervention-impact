## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# k_means.r
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


library(rasterVis)
library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(gridExtra)
library(Hmisc)

set.seed(206)

rm(list=ls())
overwrite_rotation <- F
overwrite_kmeans <- F

source("classify_functions.r")
root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")
palette <- "Paired" # color scheme for plots 

# number of singular vectors to use, from visual inspection of svd plots
input_list <- list(# tsi_rainfall = list(asia=3, americas=3) , 
                   tsi_rainfall_vector_abundance=list(africa=3)
                   )
cluster_counts <- 3:10

for (this_cov in names(input_list)){
  nvec_list <- input_list[[this_cov]]
  
  for (continent in names(nvec_list)){
    print(paste("clustering", this_cov, "for", continent))
    
    main_dir <- file.path(base_dir, continent)
    
    # rotate matrix, if needed
    rotation_fname <- file.path(main_dir, paste0("svd_rotations_", this_cov, ".csv"))
    nvecs <- nvec_list[[continent]]
    
    if (file.exists(rotation_fname) & overwrite_rotation==F){
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
    var_explained <- c()
    pdf(file.path(main_dir, paste0("k_means_", this_cov, ".pdf")), width=9, height=6)
    for (nclust in cluster_counts){
      
      cov_label <- ifelse(nchar(this_cov)==3, toupper(this_cov), capitalize(this_cov))
      
      # if k-means has already been run, just load outputs
      kmeans_dir <- file.path(main_dir, "kmeans")
      dir.create(kmeans_dir, showWarnings=F, recursive=T)
      
      k_out_fname <- file.path(kmeans_dir, paste0("k_out_", this_cov, "_", nclust, ".rdata"))
      cluster_raster_fname <- file.path(kmeans_dir, paste0("k_clusters_", this_cov, "_", nclust, ".tif"))
      random_trace_fname <- file.path(kmeans_dir, paste0("random_trace_", this_cov, "_", nclust, ".csv"))
      summary_fname <- file.path(kmeans_dir,  paste0("summary_", this_cov, "_", nclust, ".csv"))
      
      if (file.exists(k_out_fname) & overwrite_kmeans==F){
        print("k-means already run, loading outputs")
        cluster_raster <- raster(cluster_raster_fname)
        random_trace <- fread(random_trace_fname)
        summary_vals <- fread(summary_fname)
        load(k_out_fname)
        
      }else{
        print(paste("finding clusters for k of", nclust))
        k_out <- kmeans(rotation[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
        rotation[, cluster:= k_out$cluster]
        
        
        print("creating new raster")
        # load mask raster to get dimensions & extent
        temp_raster <- raster(file.path(main_dir, "rasters", "mask.tif"))

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
        all_vals[, cluster:=NULL]
        
        save(k_out, file=k_out_fname) # save last to ensure other outputs have been saved if this exists
        
      }
      
      plotlist <- NULL
      
      print("making map")
      cluster_raster <- ratify(cluster_raster)
      map_plot <- levelplot(cluster_raster, att="ID", col.regions=brewer.pal(nclust, palette),
                            xlab=NULL, ylab=NULL, scales=list(draw=F),
                            main = "", colorkey=F, margin=F)
      plotlist[[1]] <- map_plot
      
      time_series <- summary_vals[variable_name=="month"]
      time_series[, cluster:=as.factor(cluster)]
      lines <- lapply(unique(time_series$cov), function(cov_value){
        ggplot(time_series[cov==cov_value], aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
          facet_grid(cluster~.) +
          geom_line(size=2) +
          geom_line(aes(y=perc_05), size=1, linetype=2) +
          geom_line(aes(y=perc_95), size=1, linetype=2) +
          geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5) +
          # geom_line(data=random_trace[cov==cov_value], aes(group=id, color=factor(cluster)), size=0.5, alpha=0.25) + 
          scale_color_brewer(type="qual", palette=palette) +
          scale_fill_brewer(type="qual", palette=palette) + 
          scale_x_continuous(breaks=1:12, labels=1:12) +
          theme_minimal() +
          theme(legend.position="none",
                plot.title = element_text(face="bold")) +
          labs(title=paste(cov_value),
               x="Month",
               y="")
      })
      
      plotlist <- append(plotlist, lines)
      
      # if applicable, plot vector mix
      summary_species <- summary_vals[variable_name=="species"]

      if (nrow(summary_species)>0){
        vector_mix <- ggplot(summary_species, aes(x=cluster, y=mean)) +
          geom_bar(aes(fill=variable_val), stat="identity") + 
          theme_minimal() +
          theme(legend.position="bottom") +
          scale_fill_manual(values=c("#ffd92f", "#a6d854", "#8da0cb"), name="")+ 
          labs(title="Relative Species Abundance",
               x="Cluster",
               y="")
        
        vector_idx <- length(plotlist) +1 
        plotlist[[vector_idx]] <- vector_mix
        
        layout_pattern <- c(1, 1:(length(lines)+1))
        layout <- rbind(c(1, 1:(length(lines)+1)),
                        c(1, 1:(length(lines)+1)),
                        c(vector_idx, vector_idx, 2:(length(lines)+1))
                        )
        
      }else{
        layout_pattern <- c(1, 1:length(plotlist))
        layout <- rbind(layout_pattern,
                        layout_pattern)
      }
      
      print("saving")
      full_plot <- grid.arrange(grobs=plotlist, layout_matrix=layout)
      print(full_plot)
      
      # update variance explained
      var_explained <- c(var_explained, k_out$betweenss/k_out$totss)
    }
    graphics.off()
    
    # plot % variance explained across clusters
    print("plotting variance explained")
    all_var <- data.table(k=cluster_counts, var=var_explained)
    png(file=file.path(main_dir, paste0("elbow_", this_cov, ".png")))
    elbow_plot<- ggplot(all_var, aes(x=k, y=var)) +
        geom_line() +
        geom_point() +
        theme_minimal()+
        labs(x="Cluster Count",
             y="Variance Explained")
    print(elbow_plot)
    graphics.off()
  }
}

