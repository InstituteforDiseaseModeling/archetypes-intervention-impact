## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 03_cluster_plots.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# September 2018 
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
library(plotly)

set.seed(206)
theme_set(theme_minimal(base_size = 12))

rm(list=ls())

plot_vectors <- T
out_subdir <- "v5_interventions_only"

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/")

out_dir <- file.path(base_dir, "results", out_subdir)
guide <- fread(file.path(out_dir, "instructions.csv"))

# temp
guide <- guide[covariate=="itn_coverage"]

palette <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B33539", "#367A40")

for (this_continent in unique(guide$continent)){

  print(paste("plotting archetypes for", this_continent))
  
  nvecs <- unique(guide[continent==this_continent]$singular_vectors)
  
  this_in_dir <- file.path(out_dir, this_continent, "02_kmeans")
  this_out_dir <- file.path(this_in_dir, "figures")
  dir.create(this_out_dir, showWarnings=F, recursive=T)
  svd_dir <- file.path(this_in_dir, "../01_svd")

  # find cluster counts
  cluster_counts <- list.files(this_in_dir)[list.files(this_in_dir) %like% "k_out"]
  cluster_counts <- sort(as.integer(gsub("k_out_([0-9]+)_.*", "\\1", cluster_counts)))
  
  print("loading matrix rotations")
  rotation <- fread(file.path(this_in_dir, "svd_rotations.csv"))
  
  # load svd outputs to get dataset of cell values
  load(file.path(svd_dir, "svd_output.rdata")) 
  all_vals <- melt(svd_wide_datatable, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
  all_vals[, id:=as.integer(as.character(id))]
  rm(svd_out, svd_wide_datatable); gc()
  
  # full plotting
  var_explained <- c()
  pdf(file.path(this_out_dir, paste0("k_means_summary_plots.pdf")), width=8, height=4.5)
  
  print("plotting results")
  for (nclust in cluster_counts){
    
    print(paste("plotting for cluster count", nclust))
    
    if (nclust<=length(palette)){
      these_colors <- palette[1:nclust]
    }else{
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      these_colors <- gg_color_hue(nclust)
    }
    
    print("loading outputs")
    cluster_raster <- raster(file.path(this_in_dir, paste0("map_", nclust, "_cluster", ".tif")))
    random_trace <- fread(file.path(this_in_dir, paste0("random_trace_", nclust, "_cluster",  ".csv")))
    summary_vals <- fread(file.path(this_in_dir,  paste0("summary_", nclust, "_cluster", ".csv")))
    site_ids <- fread(file.path(this_in_dir,  paste0("site_ids_", nclust, "_cluster", ".csv")))
    load(file.path(this_in_dir, paste0("k_out_", nclust, "_cluster", ".rdata")))
    
    print("plotting cluster scatters")
    centers <- data.table(k_out$centers)
    setnames(centers, paste0("X", 1:nvecs), paste0("singular_vector_", 1:nvecs))
    random_trace_for_plot <- merge(unique(random_trace[, list(id, cluster)]), rotation, by="id", all.x=T)
    setnames(random_trace_for_plot, paste0("X", 1:nvecs), paste0("singular_vector_", 1:nvecs))
    
    if (nvecs==2){
      svd_cluster_plot <- ggplot(random_trace_for_plot, aes(x=singular_vector_1, y=singular_vector_2, color=as.factor(cluster))) +
        geom_point(size=3, alpha=0.8) +
        geom_point(data=centers, size=3, color="black") + 
        scale_color_manual(values=these_colors) +
        theme(legend.position = "none")
      ggsave(file.path(this_out_dir, paste0("scatter_", nclust, "_cluster", ".pdf")))
    }else if (nvecs==3){
      centers[, cluster:="centroid"]
      for_svd_plot <- rbind(random_trace_for_plot, centers, fill=T)
      for_svd_plot[, cluster:=as.factor(cluster)]
      for_svd_plot <- rbind(for_svd_plot, site_ids[, list(id, singular_vector_1, singular_vector_2, singular_vector_3, cluster="NN")])
      
      html_dir <- file.path(this_out_dir, "plotly_html")
      dir.create(html_dir, showWarnings=F, recursive=T)
      
      svd_cluster_plot <- plot_ly(for_svd_plot, x = ~singular_vector_1, y = ~singular_vector_2,
                                  color = ~cluster, colors=c(these_colors, "#000000", "#778899"),
                                  z = ~singular_vector_3, opacity=0.95) %>%
        add_markers()
      htmlwidgets::saveWidget(svd_cluster_plot, file.path(html_dir, paste0("scatter_", nclust, "_cluster", ".html")))
    }else{
      warning("Don't know what to do with this many singular vectors")
    }
    
    # plot big plot
    plotlist <- NULL
    
    # convert sites to lat-longs spatialpoints
    site_id_spoints <- xyFromCell(cluster_raster, site_ids$id, spatial=T)
    
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
      selected_sites <- all_vals[id %in% site_ids$id & cov==cov_value & variable_name=="month"]
      selected_sites <- merge(selected_sites, site_ids[, list(id, cluster)], by="id", all=T)
      selected_sites[, cluster:=as.factor(cluster)]
      
      ggplot(data, aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
        facet_grid(cluster~., scales="free_y") +
        geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5, color=NA) +
        geom_line(size=1) +
        geom_line(aes(y=perc_05), size=0.75, linetype=2) +
        geom_line(aes(y=perc_95), size=0.75, linetype=2) +
        geom_line(data=selected_sites, aes(y=cov_val), color="black", size=1) +
        scale_color_manual(values = these_colors) +
        scale_fill_manual(values = these_colors) + 
        scale_x_continuous(breaks=seq(2,12,2), labels=c("F","A","J","A","O","D"), minor_breaks=seq(1,12,2)) +
        # scale_y_continuous(breaks=ybreaks, labels=ylabs, minor_breaks=yminorbreaks, limits=ylimit) + 
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
  pdf(file=file.path(this_out_dir, "elbow_plot.pdf"))
  elbow_plot<- ggplot(all_var, aes(x=k, y=var)) +
    geom_line(size=1) +
    geom_point(shape=1, size=4) +
    theme_minimal()+
    labs(title= "", 
         x="Cluster Count",
         y="Variance Explained")
  print(elbow_plot)
  graphics.off()
}


