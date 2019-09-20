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


library(rasterVis)
library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(gridExtra)
library(latticeExtra)
library(Hmisc)
library(FNN)

set.seed(206)
theme_set(theme_minimal(base_size = 14))

rm(list=ls())
overwrite_rotation <- F
overwrite_kmeans <- T
plot_vectors <- T
covariate_type <- "no_transmission_limits"

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/")

# number of singular vectors to use, from visual inspection of svd plots
input_list <- list(# tsi_rainfall = list(asia=3, americas=3) , 
                   tsi_rainfall_relative_vector_abundance=list(africa=3)
                   # tsi_rainfall_vector_abundance_standardized=list(africa=2),
                   # tsi_rainfall_vector_abundance_rescaled=list(africa=3)
                   )

palettes <- list(
  tsi_rainfall_relative_vector_abundance= c("#F98400", "#00A08A", "#5392C2", "#902E57", "#F2AD00", "#D71B5A", "#98B548", "#8971B3", "#B03A2E", "#2ECC71"),
  tsi_rainfall_relative_vector_abundance_standardized=c("#00A08A", "#902E57", "#5392C2", "#F98400", "#F2AD00", "#D71B5A", "#98B548", "#8971B3"),
  tsi_rainfall_relative_vector_abundance_rescaled=c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B03A2E", "#2ECC71")
)


rotate_matrix <- function(nvecs, svd_dir, cov_dir, cov="tsi"){
  load(file.path(svd_dir, paste0(cov, "_svd.rdata")))
  sing_vecs <- svd_out$u[, 1:nvecs]
  
  ## multiply by original matrix to get rotations
  for_svd <- fread(file.path(cov_dir, paste0(cov, "_vals.csv")))
  print("reshaping")
  for_svd <- dcast(for_svd, cov + variable_name + variable_val ~ id)
  print("rotating")
  rotation <- data.frame(t(t(sing_vecs)%*%as.matrix(for_svd[,4:ncol(for_svd)])))
  rotation$id <- as.integer(rownames(rotation))
  rotation <- data.table(rotation)
  return(rotation)
}


cluster_counts <- 3:10

for (this_cov in names(input_list)){
  nvec_list <- input_list[[this_cov]]
  palette <- palettes[[this_cov]]
  
  for (continent in names(nvec_list)){
    print(paste("clustering", this_cov, "for", continent))
    
    cov_dir <- file.path(base_dir, "00_covariate_extraction", covariate_type, continent)
    svd_dir <- file.path(base_dir, "01_svd", continent)
    main_dir <- file.path(base_dir, "02_kmeans", continent, this_cov)
    dir.create(main_dir, showWarnings=F, recursive=T)
    
    # rotate matrix, if needed
    rotation_fname <- file.path(svd_dir, paste0("svd_rotations_", this_cov, ".csv"))
    nvecs <- nvec_list[[continent]]
    
    if (file.exists(rotation_fname) & overwrite_rotation==F){
      print("loading matrix rotations")
      rotation <- fread(rotation_fname)
    }else{
      print("finding matrix rotations")
      rotation <- rotate_matrix(nvecs, svd_dir, cov_dir, this_cov)
      write.csv(rotation, rotation_fname, row.names=F)
    }
    
    all_vals <- fread(file.path(cov_dir, paste0(this_cov, "_vals.csv")))
    setnames(all_vals, "value", "cov_val")
    
    # k-means and plotting
    var_explained <- c()
    pdf(file.path(main_dir, paste0("k_means_", this_cov, ".pdf")), width=8, height=4.5)
    for (nclust in cluster_counts){
      
      cov_label <- ifelse(nchar(this_cov)==3, toupper(this_cov), capitalize(this_cov))
      
      # if k-means has already been run, just load outputs
      k_out_fname <- file.path(main_dir, paste0("k_out_", this_cov, "_", nclust, ".rdata"))
      cluster_raster_fname <- file.path(main_dir, paste0("k_clusters_", this_cov, "_", nclust, ".tif"))
      random_trace_fname <- file.path(main_dir, paste0("random_trace_", this_cov, "_", nclust, ".csv"))
      summary_fname <- file.path(main_dir,  paste0("summary_", this_cov, "_", nclust, ".csv"))
      site_id_fname <- file.path(main_dir,  paste0("site_ids", this_cov, "_", nclust, ".csv"))
      
      if (file.exists(k_out_fname) & overwrite_kmeans==F){
        print("k-means already run, loading outputs")
        cluster_raster <- raster(cluster_raster_fname)
        random_trace <- fread(random_trace_fname)
        summary_vals <- fread(summary_fname)
        site_ids <- fread(site_id_fname)
        load(k_out_fname)
        
      }else{
        print(paste("finding clusters for k of", nclust))
        k_out <- kmeans(rotation[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
        rotation[, cluster:= k_out$cluster]
        
        
        print("creating new raster")
        # load mask raster to get dimensions & extent
        temp_raster <- raster(file.path(cov_dir, "rasters", "mask.tif"))

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
        
        print("finding centroids for representative sites")

        # bit of matrix transposition to appropriately map a raster ID to a matrix cell
        for_centroids <- rep(NA, ncell(temp_raster))
        for_centroids[rotation$id] <- rotation$id
        for_centroids <- matrix(for_centroids, nrow=nrow(temp_raster))
        for_centroids <- base::t(for_centroids)
        cell_id_map <- data.table(new_id = which(!is.na(for_centroids)),
                              id = for_centroids[!is.na(for_centroids)])
        
        # find cell IDs of centroids and map them to the appropriate locations
        site_ids <- rotation[get.knnx(rotation[,1:nvecs], k_out$centers, k=1)$nn.index]
        site_ids <- merge(site_ids, cell_id_map, by="id", all.x=T)
        write.csv(site_ids, site_id_fname, row.names=F)
        
        save(k_out, file=k_out_fname) # save last to ensure other outputs have been saved if this exists
        
      }
      
      plotlist <- NULL
      these_colors <- palette[1:nclust]
      
      # convert sites to lat-longs spatialpoints
      site_id_spoints <- xyFromCell(cluster_raster, site_ids$new_id, spatial=T)
      
      print("making map")
      cluster_raster <- ratify(cluster_raster)
      map_plot <- levelplot(cluster_raster, att="ID", col.regions=these_colors,
                            xlab=NULL, ylab=NULL, scales=list(draw=F),
                            main = "", colorkey=F, margin=F)   +
        latticeExtra::layer(sp.points(site_id_spoints), theme = simpleTheme(col = "black",
                                                                            cex=2))
      plotlist[[1]] <- map_plot
      
      time_series <- summary_vals[variable_name=="month"]
      time_series[, cluster:=as.factor(cluster)]
      lines <- lapply(unique(time_series$cov), function(cov_value){
        
        data <- time_series[cov==cov_value]
        
        if(max(data$perc_95)>1){
          ybreaks <- c(0, 600)
          ylabs <- c("0", "600")
          yminorbreaks <- c(100, 200, 300, 400, 500)
          ylimit <- NULL
        }else if(max(data$perc_95)<0.75){
          ybreaks <- c(0, 0.75)
          ylabs <- c("0", "0.75")
          yminorbreaks <- c(0.25, 0.5)
          ylimit <- c(0, 0.75)
        }else{
          ybreaks <- c(0, 0.75)
          ylabs <- c("0", "0.75")
          yminorbreaks <- c(0.25, 0.5)
          ylimit <- NULL
        }
        
        ggplot(data, aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
          facet_grid(cluster~.) +
          geom_line(size=1) +
          geom_line(aes(y=perc_05), size=0.75, linetype=2) +
          geom_line(aes(y=perc_95), size=0.75, linetype=2) +
          geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5, color=NA) +
          scale_color_manual(values = these_colors) +
          scale_fill_manual(values = these_colors) + 
          scale_x_continuous(breaks=seq(2,12,2), labels=c("F","A","J","A","O","D"), minor_breaks=seq(1,12,2)) +
          scale_y_continuous(breaks=ybreaks, labels=ylabs, minor_breaks=yminorbreaks, limits=ylimit) + 
          theme(legend.position="none",
                plot.title = element_text(size=16),
                strip.background = element_blank(),
                strip.text.y = element_blank()) +
          labs(title=ifelse(nchar(cov_value)==3, toupper(cov_value), capitalize(cov_value)),
               x="",
               y="")
      })
      
      plotlist <- append(plotlist, lines)
      
      # if applicable, plot vector mix
      summary_species <- summary_vals[variable_name=="species"]

      if (nrow(summary_species)>0 & plot_vectors==T){
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
    pdf(file=file.path(main_dir, paste0("elbow_", this_cov, ".pdf")))
    elbow_plot<- ggplot(all_var, aes(x=k, y=var)) +
        geom_line(size=1) +
        geom_point(shape=1, size=4) +
        theme_minimal()+
        labs(title= "", # paste("Elbow Plot,", this_cov),
             x="Cluster Count",
             y="Variance Explained")
    print(elbow_plot)
    graphics.off()
  }
}

