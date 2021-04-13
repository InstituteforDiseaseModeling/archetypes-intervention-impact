## -----------------------------------------------------------------------------------------------------------------
# arch_paper_plots.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# April 2021 
# 
# Plots for archetypes/intervention impact paper, draft here:  
# https://www.overleaf.com/read/hqbdxggcwvwv
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(rasterVis)
library(data.table)
library(ggplot2)
library(gridExtra)
library(latticeExtra)

rm(list = ls())

palette <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B33539", "#367A40")

### Setup -----------------------------------------------------------------------------------------------------------------------

main_dir <- "~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/"
arch_dir <- file.path(main_dir, "archetypes/results/v4_era5_bounded_transmission/africa")
ii_dir <- file.path(main_dir, "intervention_impact/20191009_mega_era5_new_arch")
sensitivity_dir <- file.path(main_dir, "intervention_impact/20191218_site_sensitivity")

africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"
source("~/repos/malaria-atlas-project/intervention_impact/visualize_results/map_ii_functions.r")

out_dir <- file.path(main_dir, "writing_and_presentations/ii_paper/figures/raw_figs")


### Load and Prep Data -----------------------------------------------------------------------------------------------------------------------

get_smooth <- function(x, y){
  if (max(y)<0.05){
    return(y)
  }else{
    lo <- loess(y[y>0]~x[y>0])
    predictions <- c(y[y==0], predict(lo))
    return(pmax(predictions, rep(0, length(predictions))))
  }
}

smooth_data <- function(input){
  minmaxes_smooth <- lapply(unique(input$Site_Name), function(site_name){
    sub_list <- lapply(unique(input$int_id), function(int_name){
      subset <- input[Site_Name==site_name & int_id==int_name]
      subset[, smooth_min:= get_smooth(mean_initial, min_final)]
      subset[, smooth_max:= get_smooth(mean_initial, max_final)]
      subset[, smooth_mean:= get_smooth(mean_initial, mean_final)]
    })
    sub_list <- rbindlist(sub_list)
  })
  minmaxes_smooth <- rbindlist(minmaxes_smooth)
  return(minmaxes_smooth)
}


ii_in_dir <- file.path(ii_dir, "results/clean")
full_impact <- fread(file.path(ii_in_dir, "full_impact.csv"))
summary_impact <- fread(file.path(ii_in_dir, "summary_impact.csv"))
smooth_impact <- smooth_data(summary_impact)
smooth_impact[, Site_Name:= factor(Site_Name)]

impact_brick <- brick(file.path(ii_dir, "results/rasters/pfpr_Africa.tif"))


### Shapefile  -----------------------------------------------------
africa_shp <- readOGR(africa_shp_dir)
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)


### Explainer Figure -----------------------------------------------------------------------------------------------------------------------

## Part 1: 4-Site Map
cluster_counts <- c(4, 10)
for (nclust in cluster_counts){
  print(paste("generating archetypes plot for cluster count", nclust))
  
  # load k-means results
  cluster_in_dir <- file.path(arch_dir, "02_kmeans")
  cluster_raster <- raster(file.path(cluster_in_dir, paste0("map_", nclust, "_cluster", ".tif")))
  random_trace <- fread(file.path(cluster_in_dir, paste0("random_trace_", nclust, "_cluster",  ".csv")))
  summary_vals <- fread(file.path(cluster_in_dir,  paste0("summary_", nclust, "_cluster", ".csv")))
  site_ids <- fread(file.path(cluster_in_dir,  paste0("site_ids_", nclust, "_cluster", ".csv")))
  
  # load svd outputs to get dataset of cell values
  load(file.path(cluster_in_dir, "../01_svd/svd_output.rdata")) 
  all_vals <- melt(svd_wide_datatable, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
  all_vals[, id:=as.integer(as.character(id))]
  rm(svd_wide_datatable); gc()
  
  # rename covariates
  summary_vals[, cov:=factor(cov,
                             levels=c("air_temp_era5", "precip_era5", "relative_vector_abundance"),
                             labels=c("Air Temperature", "Precipitation", "Rel. Vector\nAbundance"))]
  all_vals[, cov:=factor(cov,
                         levels=c("air_temp_era5", "precip_era5", "relative_vector_abundance"),
                         labels=c("Air Temperature", "Precipitation", "Rel. Vector\nAbundance"))]
  
  
  # plot big plot
  plotlist <- NULL
  these_colors <- palette[1:nclust]
  
  # convert sites to lat-longs spatialpoints
  site_id_spoints <- xyFromCell(cluster_raster, site_ids$id, spatial=T)
  
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
      theme_minimal() +
      theme(legend.position="none",
            plot.title = element_text(size=16),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      labs(title=cov_value,
           x="",
           y="")
  })
  
  plotlist <- append(plotlist, lines)
  
  
  #plot vector mix
  vector_props <- all_vals[variable_name=="species" & id %in% site_ids$id]
  vector_props <- merge(vector_props, site_ids[, list(id, cluster)], by="id", all=T)
  vector_props[, cluster:=as.factor(cluster)]

  vector_props[, fraction:=cov_val]
  vector_props[, ymax:=cumsum(fraction), by="cluster"] # Compute the cumulative percentages (top of each rectangle)
  vector_props[, ymin:=ymax-fraction] # Compute the bottom of each rectangle
  
  vector_mix <- ggplot(vector_props, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(aes(fill=variable_val)) + 
    theme_void() +
    theme(legend.position="bottom",
          strip.background = element_blank(),
          strip.text = element_blank()) +
    scale_fill_manual(values=c("#ffd92f", "#a6d854", "#8da0cb"), name="")+ 
    facet_grid(cluster ~ .) + 
    coord_polar(theta = "y") +
    xlim(c(2,4)) +
    labs(title="Vector Abundance",
         x="Cluster",
         y="")
  
  vector_idx <- length(plotlist) +1 
  plotlist[[vector_idx]] <- vector_mix
  
  layout <- matrix(c(1, 1:(length(lines)+1), vector_idx), nrow=1)
  
  pdf(file.path(out_dir, paste0("archs_", nclust,".pdf")), width=8, height=4.5)
    grid.arrange(grobs=plotlist, layout_matrix=layout)
  graphics.off()
  
}


## Part 2: Intervention Impact Curves
int_to_use <- 79 # 79 is irs: 0%, itn: 60%, al_cm: 60% 

pdf(file.path(out_dir, "example_curves.pdf"), width=8, height=4.5)
ggplot(smooth_impact[int_id==int_to_use & Site_Name %in% c(6, 5, 8, 10)], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.5)+
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=Site_Name, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=Site_Name, group=Site_Name), size=1.25) +
  scale_color_manual(values=these_colors, name="Site") +
  scale_fill_manual(values=these_colors, name="Site") +
  xlim(0,0.8) +
  ylim(0,0.8) +
  facet_grid(. ~ Site_Name) + 
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  coord_fixed() +
  labs(x="Initial Prevalence",
       y="Final Prevalence",
       title="")
graphics.off()


## Part 3: Example output maps 
color_vals <- generate_full_pal()
pdf(file.path(out_dir, "example_maps.pdf"), width=8, height=4.5)
# first layer is PfPR in 2000, other layers are interventions in order
print(levelplot(impact_brick[[c(1, int_to_use+1)]], par.settings=rasterTheme(color_vals$pal), at=color_vals$breaks,
                xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F)) +
        latticeExtra::layer(sp.polygons(africa_shp)))
graphics.off()



### 10-Site Archetype Figure -----------------------------------------------------------------------------------------------------------------------

## see explainer figure code

### Intervention Impact by Biting Intensity for Vector and Non-vector Interventions  -----------------------------------------------------------------------------------------------------------------------
vector_data <- fread(file.path(ii_dir, "input/vector/vector_proportions.csv"))
anthro_endo <- fread(file.path(ii_dir, "results/tables_for_paper/vector_anthro_endo.csv"))
vector_data <- melt(vector_data, id.vars=c("id", "continent", "lat", "lon"), variable.name = "species", value.name="fraction")
vector_data <- merge(vector_data, anthro_endo)

site_anthro_endo <- vector_data[, list(anthro=sum(anthro*fraction), endo=sum(endo*fraction)), by=list(id, continent)]
site_anthro_endo[, human_indoor:= round(anthro*endo, 2)]
setnames(site_anthro_endo, "id", "Site_Name")
smooth_impact <- merge(smooth_impact, site_anthro_endo, all.x=T)
smooth_impact[, human_indoor:=factor(human_indoor)]

these_colors <- c(palette, "black")

pdf(file.path(out_dir, "biting_intensity.pdf"), width=11, height=7)
ggplot(smooth_impact[!Site_Name %in% c(1, 11, 12) & int_id %in% c(101, 25)], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.5)+
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=human_indoor, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=human_indoor, group=Site_Name), size=1.25) +
  scale_color_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  scale_fill_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  xlim(0,0.85) +
  ylim(0,0.85) +
  facet_grid(.~label) +
  coord_fixed() +
  labs(x="Initial Prevalence",
       y="Final Prevalence",
       title="")
graphics.off()




### Example Intervention Packages and Maps -----------------------------------------------------------------------------------------------------------------------


### Sensistivity Analysis Figure -----------------------------------------------------------------------------------------------------------------------
sensitivity_in_dir <- file.path(sensitivity_dir, "results/clean")
#sensitivity_full_impact <- fread(file.path(sensitivity_in_dir, "full_impact.csv"))
sensitivity_summary_impact <- fread(file.path(sensitivity_in_dir, "summary_impact.csv"))
sensitivity_smooth <- smooth_data(sensitivity_summary_impact)


### Supplement: Covariate Normalization -----------------------------------------------------------------------------------------------------------------------


### Supplement: SVD Vectors -----------------------------------------------------------------------------------------------------------------------


### Supplement: K-means elbow -----------------------------------------------------------------------------------------------------------------------


### Supplement: All cluster maps -----------------------------------------------------------------------------------------------------------------------


### Supplement: All package plots and maps -----------------------------------------------------------------------------------------------------------------------