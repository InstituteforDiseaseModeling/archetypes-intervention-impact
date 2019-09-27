## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 03_cluster_plots.r
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
library(plotly)

set.seed(206)
theme_set(theme_minimal(base_size = 14))

rm(list=ls())

plot_vectors <- T
out_subdir <- "v1_original_megatrends"
root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/")
out_dir <- file.path(base_dir, "results", out_subdir)
guide <- fread(file.path(out_dir, "instructions.csv"))

palette <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B03A2E")

this_continent <- "africa"
nclust <- 6
these_colors <- palette[1:nclust]
nvecs <- 3

this_in_dir <- file.path(out_dir, this_continent, "02_kmeans")
this_out_dir <- file.path(this_in_dir, "figures")
svd_dir <- file.path(this_in_dir, "../01_svd")

print("loading matrix rotations")
rotation <- fread(file.path(this_in_dir, "svd_rotations.csv"))

# load svd outputs to get dataset of cell values
load(file.path(svd_dir, "svd_output.rdata")) 
all_vals <- melt(svd_wide_datatable, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
all_vals[, id:=as.integer(as.character(id))]
rm(svd_out, svd_wide_datatable); gc()

# full plotting

print("plotting results")

print("loading outputs")
cluster_raster <- raster(file.path(this_in_dir, paste0("map_", nclust, "_cluster", ".tif")))
random_trace <- fread(file.path(this_in_dir, paste0("random_trace_", nclust, "_cluster",  ".csv")))
summary_vals <- fread(file.path(this_in_dir,  paste0("summary_", nclust, "_cluster", ".csv")))
load(file.path(this_in_dir, paste0("k_out_", nclust, "_cluster", ".rdata")))

# NEW: load custom site locations from ad-hoc selection
custom_site_ids <- fread("~/repos/malaria-atlas-project/intervention_impact/run_simulations/input_files/site_details.csv")
custom_site_ids <- custom_site_ids[continent=="Africa", list(name, lat, lon, country)]
custom_site_ids$id <- cellFromXY(cluster_raster, as.matrix(custom_site_ids[, list(lon, lat)]))

# MANUAL map to id vals
cluster_map <- data.table(name=c("djibo", "gode", "aba", "kananga", "kasama", "moine"),
                          cluster=c(4,2,1,3,6,5))
custom_site_ids <- merge(custom_site_ids, cluster_map, by="name", all=T)

print("plotting cluster scatters")
centers <- data.table(k_out$centers)
centers[, cluster:="centroid"]
random_trace_for_plot <- merge(unique(random_trace[, list(id, cluster)]), rotation, by="id", all.x=T)
custom_site_rotation <- rotation[id %in% custom_site_ids$id]
custom_site_rotation[, cluster:="custom ID"]

for_svd_plot <- rbind(random_trace_for_plot, centers, fill=T)
for_svd_plot <- rbind(for_svd_plot, custom_site_rotation)
for_svd_plot[, cluster:=as.factor(cluster)]


html_dir <- file.path(this_out_dir, "plotly_html")

svd_cluster_plot <- plot_ly(for_svd_plot, x = ~X1, y = ~X2, color = ~cluster, colors=c(these_colors, "#000000", "#2F4F4F"), z = ~X3, opacity=0.85) %>%
  add_markers()
htmlwidgets::saveWidget(svd_cluster_plot, file.path(html_dir, paste0("CUSTOM_scatter_", nclust, "_cluster", ".html")))


# plot big plot
plotlist <- NULL

pdf(file.path(this_out_dir, paste0("CUSTOM_k_means_6cluster.pdf")), width=8, height=4.5)
# convert sites to lat-longs spatialpoints
site_id_spoints <- xyFromCell(cluster_raster, custom_site_ids$id, spatial=T) 
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
  selected_sites <- all_vals[id %in% custom_site_ids$id & cov==cov_value & variable_name=="month"]
  selected_sites <- merge(selected_sites, custom_site_ids[, list(id, cluster)], by="id", all=T)
  selected_sites[, cluster:=as.factor(cluster)]
  
  ggplot(data, aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
    facet_grid(cluster~.) +
    geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5, color=NA) +
    geom_line(size=1) +
    geom_line(aes(y=perc_05), size=0.75, linetype=2) +
    geom_line(aes(y=perc_95), size=0.75, linetype=2) +
    geom_line(data=selected_sites, aes(y=cov_val), color="black", size=1) +
    scale_color_manual(values = these_colors) +
    scale_fill_manual(values = these_colors) + 
    scale_x_continuous(breaks=seq(2,12,2), labels=c("F","A","J","A","O","D"), minor_breaks=seq(1,12,2)) +
    theme(legend.position="none",
          plot.title = element_text(size=16),
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    labs(title=ifelse(nchar(cov_value)==3, toupper(cov_value), capitalize(cov_value)),
         x="",
         y="")
})

plotlist <- append(plotlist, lines)

# plot vector mix
summary_species <- summary_vals[variable_name=="species"]
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


print("saving")
full_plot <- grid.arrange(grobs=plotlist, layout_matrix=layout)
print(full_plot)

graphics.off()




