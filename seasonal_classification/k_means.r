library(rasterVis)
library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(gridExtra)
library(Hmisc)

rm(list=ls())

source("classify_functions.r")
base_dir <- file.path(Sys.getenv("USERPROFILE"), "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")
palette <- "Paired"

# number of singular vectors to use, from visual inspection of svd plots
nvec_list <- list(africa=2, asia=3, americas=2)
# nvec_list <- list(africa_orig=2)
cov_list <- c("tsi")

for (this_cov in cov_list){
  for (continent in names(nvec_list)){
    print(paste("clustering", this_cov, "for", continent))
    
    main_dir <- file.path(base_dir, continent)
    
    # rotate matrix, if needed
    rotation_fname <- file.path(main_dir, paste0("svd_rotations_", this_cov, ".csv"))
    
    if (file.exists(rotation_fname)){
      rotation <- fread(rotation_fname)
    }else{
      print("finding matrix rotations")
      nvecs <- nvec_list[[continent]]
      rotation <- rotate_matrix(nvecs, main_dir, this_cov)
      write.csv(rotation, rotation_fname, row.names=F)
    }
    
    all_vals <- fread(file.path(main_dir, paste0(this_cov, "_vals.csv")))
    setnames(all_vals, this_cov, "cov_val")
    
    # Plotting
    pdf(file.path(main_dir, paste0("k_means_", this_cov, ".pdf")), width=9, height=6)
    for (nclust in 3:10){
      print(paste("finding clusters for k of", nclust))
      k_out <- kmeans(rotation[, 1:2], centers=nclust, algorithm = "MacQueen", iter.max=100)
      rotation[, cluster:= k_out$cluster]
      
      cov_label <- ifelse(nchar(this_cov)==3, toupper(this_cov), capitalize(this_cov))
      
      print("making map")
      new <- matrix(nrow=1609, ncol=1945)
      new[rotation$id] <- rotation$cluster
      new <- ratify(raster(new))
      map_plot <- levelplot(new, att="ID", col.regions=brewer.pal(nclust, palette),
                            xlab=NULL, ylab=NULL, scales=list(draw=F),
                            main = paste("Mapped", cov_label, "Clusters"), colorkey=F)
      
      print("finding mean time series")
      all_vals <- merge(all_vals, rotation[, list(id, cluster)], by="id", all=T)
      summary_vals <- all_vals[, list(cov_val=mean(cov_val)), by=list(cluster, month)]
      all_vals[, cluster:=NULL]
      
      lines <- ggplot(summary_vals, aes(x=month, y=cov_val)) +
        geom_line(aes(color=factor(cluster)), size=2) +
        facet_grid(cluster~.) +
        scale_color_brewer(type="qual", palette=palette) +
        scale_x_continuous(breaks=1:12, labels=1:12) +
        theme_minimal() +
        theme(legend.position="none",
              plot.title = element_text(face="bold")) +
        labs(title=paste("Mean", cov_label, "Trend"),
             x="Month",
             y=cov_label)
      
      layout <- rbind(c(1,1,2),
                      c(1,1,2))
      
      print("saving")
      full_plot <- grid.arrange(map_plot, lines, layout_matrix=layout)
      print(full_plot)
      
    }
    graphics.off()
    
  }
}




