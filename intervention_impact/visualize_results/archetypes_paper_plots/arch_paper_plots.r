## -----------------------------------------------------------------------------------------------------------------
# arch_paper_plots.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# April 2021 
# 
# Plots for archetypes/intervention impact paper, draft here:  
# https://www.overleaf.com/read/hqbdxggcwvwv
## -----------------------------------------------------------------------------------------------------------------------

library(rgdal)
library(stringr)
library(raster)
library(rasterVis)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(latticeExtra)
library(bit64)

rm(list = ls())

palette <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B33539", "#367A40", "#063970", "#EACB44", "#8FBC94", "#75CCD0")

### Setup -----------------------------------------------------------------------------------------------------------------------

main_dir <- "~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/"

results_type <- "arch"
plot_sens <- T

if (results_type=="arch"){
  arch_dir <- file.path(main_dir, "archetypes/results/v4_era5_bounded_transmission/africa")
  ii_dir <- file.path(main_dir, "intervention_impact/20191009_mega_era5_new_arch")
  impact_raster_fname <- file.path(ii_dir, "results/rasters/pfpr_Africa.tif")
  final_nclust <- 10 # number of clusters in actual analysis
  sites_to_use <- 2:10
  orig_cov_names <- c("air_temp_era5", "precip_era5", "relative_vector_abundance")
  pretty_cov_names <- c("Air Temperature", "Precipitation", "Rel. Vector\nAbundance")
  out_dir <- file.path(main_dir, "writing_and_presentations/ii_paper/figures/raw_figs")
}else if (results_type=="mega"){
  plot_sens <- F
  arch_dir <- file.path(main_dir, "archetypes/results/v1_original_megatrends/africa")
  ii_dir <- file.path(main_dir, "intervention_impact/20191008_replicate_megatrends")
  impact_raster_fname <- file.path(ii_dir, "results/rasters/pfpr_Africa.tif")
  final_nclust <- 6 # no clusters in actual analysis
  sites_to_use <- 1:6
  orig_cov_names <- c("tsi", "rainfall_chirps", "relative_vector_abundance")
  pretty_cov_names <- c("TSI", "Rainfall", "Rel. Vector\nAbundance")
  out_dir <- file.path(main_dir, "writing_and_presentations/megatrends/malj_paper/figures/raw_figs")
}
dir.create(out_dir, showWarnings = F)

sensitivity_dir <- file.path(main_dir, "intervention_impact/20191218_site_sensitivity")

africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"
cov_dir <- file.path(main_dir, "archetypes/covariates/with_transmission_limits/africa")
source("~/repos/malaria-atlas-project/intervention_impact/visualize_results/map_ii_functions.r")



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
summary_impact[, Site_Name:= as.character(Site_Name)]
smooth_impact <- smooth_data(summary_impact)
smooth_impact[, Site_Name:= factor(Site_Name)]

smooth_impact[, label:= factor(paste0(str_pad(int_id, 3, side="left", pad="0"), ": ", gsub(" ", "", label)))]

impact_brick <- brick(impact_raster_fname)
names(impact_brick) <- c("baseline", as.character(unique(smooth_impact$label)))

### Shapefile  -----------------------------------------------------
africa_shp <- readOGR(africa_shp_dir)
africa_dt <- data.table(fortify(africa_shp, region = "COUNTRY_ID"))
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)


### 10-site color map
model_archs_colormap <- fread(file.path(ii_dir, "results", "clean", "cluster_color_map.csv"))
model_archs_colormap <- model_archs_colormap[order(ns_order)]
if ("Site_Name" %in% names(model_archs_colormap)){
  model_archs_colormap[, Site_Name:=as.factor(Site_Name)]
}

ns_palette <- model_archs_colormap$color


### Archetype Figures -----------------------------------------------------------------------------------------------------------------------


raster_plotlist <- NULL
raster_plotlist_idx <- 1

for (nclust in 3:14){
  print(paste("generating archetypes plot for cluster count", nclust))
  
  # load k-means results
  cluster_in_dir <- file.path(arch_dir, "02_kmeans")
  cluster_raster <- raster(file.path(cluster_in_dir, paste0("map_", nclust, "_cluster", ".tif")))
  random_trace <- fread(file.path(cluster_in_dir, paste0("random_trace_", nclust, "_cluster",  ".csv")))
  summary_vals <- fread(file.path(cluster_in_dir,  paste0("summary_", nclust, "_cluster", ".csv")))
  site_ids <- fread(file.path(cluster_in_dir,  paste0("site_ids_", nclust, "_cluster", ".csv")))
  site_ids <- site_ids[order(-latitude)]
  site_ids[, ns_order:= as.integer(row.names(site_ids))]
  
  custom_site_id_fname <- file.path(cluster_in_dir,  paste0("site_ids_", nclust, "_cluster_CUSTOM", ".csv"))
  if (file.exists(custom_site_id_fname)){
    custom_sites <- T
    custom_site_ids <- fread(custom_site_id_fname)
    custom_site_ids <- merge(custom_site_ids, model_archs_colormap[, list(cluster, ns_order)], by="cluster", all.x=T)
    # custom_site_ids <- custom_site_ids[order(-latitude)]
    # custom_site_ids[, ns_order:= as.integer(row.names(custom_site_ids))]
    custom_site_ids$id <- cellFromXY(cluster_raster, as.matrix(custom_site_ids[, list(longitude, latitude)]))
  }else{
    custom_sites <- F
  }
  
  # load svd outputs to get dataset of cell values
  load(file.path(cluster_in_dir, "../01_svd/svd_output.rdata")) 
  all_vals <- melt(svd_wide_datatable, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
  all_vals[, id:=as.integer(as.character(id))]
  rm(svd_wide_datatable); gc()
  
  # rename covariates
  summary_vals[, cov:=factor(cov,
                             levels=orig_cov_names,
                             labels=pretty_cov_names)]
  all_vals[, cov:=factor(cov,
                         levels=orig_cov_names,
                         labels=pretty_cov_names)]
  
  
  # plot big plot
  plotlist <- NULL
  these_colors <- palette[1:nclust]
  
  # convert sites to lat-longs spatialpoints
  site_id_spoints <- xyFromCell(cluster_raster, site_ids$id, spatial=T)
  if (custom_sites){
    custom_site_id_spoints <- xyFromCell(cluster_raster, custom_site_ids$id, spatial=T)
  }
  
  cluster_dt <- data.table(rasterToPoints(cluster_raster))
  names(cluster_dt) <- c("long", "lat", "value")
  cluster_dt[, value:= as.factor(value)]
  
  cluster_plot <- ggplot() +
    geom_raster(data = cluster_dt, aes(fill = value, y = lat, x = long)) +
    geom_path(data = africa_dt, aes(x = long, y = lat, group = group), color = "black", size = 0.3) + 
    scale_fill_manual(values= these_colors) +
    coord_equal(xlim = c(-18, 52), ylim = c(-35, 38)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in"), legend.position = "none")
  
  # cluster_raster <- ratify(cluster_raster)
  # map_plot <- levelplot(cluster_raster, att="ID", col.regions=these_colors,
  #                       xlab=NULL, ylab=NULL, scales=list(draw=F),
  #                       main = "", colorkey=F, margin=F)   +
  #   latticeExtra::layer(sp.points(site_id_spoints), theme = simpleTheme(col = "black",
  #                                                                       cex=2))
  
  if (nclust == final_nclust){
    cluster_plot_for_key <- cluster_plot
  }
  
  cluster_plot <- cluster_plot + 
    geom_point(data=data.table(site_id_spoints@coords), aes(x=x, y=y), color="black", shape=3, size=3)
  
  if (custom_sites){
    cluster_plot <- cluster_plot + 
      geom_point(data=data.table(custom_site_id_spoints@coords), aes(x=x, y=y), color="black", shape=1, size=5)
  }
  
  plotlist[[1]] <- cluster_plot
  raster_plotlist[[raster_plotlist_idx]] <- cluster_plot
  
  time_series <- summary_vals[variable_name=="month"]
  time_series  <- merge(time_series, site_ids[, list(cluster, ns_order)], by="cluster", all.x=T)
  time_series[, cluster:=as.factor(cluster)]
  
  
  lines <- lapply(unique(time_series$cov), function(cov_value){
    
    data <- time_series[cov==cov_value]
    selected_sites <- all_vals[id %in% site_ids$id & cov==cov_value & variable_name=="month"]
    selected_sites <- merge(selected_sites, site_ids[, list(id, cluster, ns_order)], by="id", all=T)
    selected_sites[, cluster:=as.factor(cluster)]
    
    if(custom_sites){
      custom_selected_sites <- all_vals[id %in% custom_site_ids$id & cov==cov_value & variable_name=="month"]
      custom_selected_sites <- merge(custom_selected_sites, custom_site_ids[, list(id, cluster, ns_order)], by="id", all=T)
      custom_selected_sites[, cluster:=as.factor(cluster)]
      
      selected_sites[, ns_order:=NULL]
      selected_sites <- merge(selected_sites, unique(custom_selected_sites[, list(cluster, ns_order)]), by="cluster", all=T)
      
      data[, ns_order:=NULL]
      data <- merge(data, unique(custom_selected_sites[, list(cluster, ns_order)]), by="cluster", all=T)
      
    }
    
    if (max(data$perc_95)>1){
      y_breaks <- seq(0, round(max(data$perc_95)), length.out=3)
      y_labs <- c("0", as.character(round(max(data$perc_95)/2)), as.character(round(max(data$perc_95))))
      y_lims <- c(0, round(max(data$perc_95)))
    }else{
      y_breaks <- c(0,0.5,1)
      y_labs <- c("0", "0.5", "1")
      y_lims <- c(0,1)
    }

    
    line_plot <- ggplot(data, aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
                  facet_grid(ns_order~.) +
                  geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5, color=NA) +
                  geom_line(size=0.5) +
                  geom_line(aes(y=perc_05), size=0.25, linetype=2) +
                  geom_line(aes(y=perc_95), size=0.25, linetype=2) +
                  geom_line(data=selected_sites, aes(y=cov_val), color="black", size=0.5) +
                  scale_color_manual(values = these_colors) +
                  scale_fill_manual(values = these_colors) + 
                  scale_x_continuous(breaks=seq(2,12,2), labels=c("F","A","J","A","O","D"), minor_breaks=seq(1,12,2)) +
                  scale_y_continuous(breaks=y_breaks, labels=y_labs, limits=y_lims) +
                  theme_minimal(base_size=8) +
                  theme(legend.position="none",
                        # plot.title = element_text(size=8),
                        strip.background = element_blank(),
                        strip.text.y = element_blank()) +
                  labs(title=cov_value,
                       x="",
                       y="")
    
    if (custom_sites){
      line_plot <- line_plot + 
        geom_line(data=custom_selected_sites, aes(y=cov_val),  color="black", size=0.75, linetype="dashed") #
    }
    
    return(line_plot)
    
  })
  
  plotlist <- append(plotlist, lines)
  
  
  #plot vector mix
  
  if (custom_sites){
    vector_props <- all_vals[variable_name=="species" & id %in% custom_site_ids$id]
    vector_props <- merge(vector_props, custom_site_ids[, list(id, cluster, ns_order)], by="id", all=T)
    vector_props[, cluster:=as.factor(cluster)]
  }else{
    vector_props <- all_vals[variable_name=="species" & id %in% site_ids$id]
    vector_props <- merge(vector_props, site_ids[, list(id, cluster, ns_order)], by="id", all=T)
    vector_props[, cluster:=as.factor(cluster)]
  }
  
  vector_props[, fraction:=cov_val]
  vector_props[, ymax:=cumsum(fraction), by="cluster"] # Compute the cumulative percentages (top of each rectangle)
  vector_props[, ymin:=ymax-fraction] # Compute the bottom of each rectangle
  vector_props[, variable_val:= factor(variable_val, levels=c("arabiensis", "coluzzii_gambiae", "funestus_subgroup"), labels=c("ara", "gam", "fun"))]
  
  vector_mix <- ggplot(vector_props, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(aes(fill=variable_val)) + 
    theme_void(base_size=8) +
    theme(legend.position="bottom",
          strip.background = element_blank(),
          # plot.title = element_text(size=8),
          strip.text = element_blank()) +
    scale_fill_manual(values=c("#ffd92f", "#a6d854", "#8da0cb"), name="")+ 
    facet_grid(ns_order ~ .) + 
    coord_polar(theta = "y") +
    xlim(c(2,4)) +
    labs(title="Species",
         x="Cluster",
         y="")
  
  vector_idx <- length(plotlist) +1 
  plotlist[[vector_idx]] <- vector_mix
  
  layout <- matrix(c(1, 1:(length(lines)+1), vector_idx), nrow=1)
  
  pdf(file.path(out_dir, paste0("archs_", nclust,".pdf")), width=8, height=4.5)
  grid.arrange(grobs=plotlist, layout_matrix=layout)
  graphics.off()
  
  raster_plotlist_idx <- raster_plotlist_idx +1
}


pdf(file.path(out_dir, paste0("arch_maps_all.pdf")), width=7, height=9)
  do.call("grid.arrange", c(raster_plotlist, ncol=3))
graphics.off()


### Intervention Impact by Biting Intensity for Vector and Non-vector Interventions  -----------------------------------------------------------------------------------------------------------------------
vector_data <- fread(file.path(ii_dir, "input/vector/vector_proportions.csv"))
anthro_endo <- fread(file.path(ii_dir, "results/tables_for_paper/vector_anthro_endo.csv"))
vector_data <- melt(vector_data, id.vars=c("id", "continent", "lat", "lon"), variable.name = "species", value.name="fraction")
vector_data <- merge(vector_data, anthro_endo)

site_anthro_endo <- vector_data[, list(anthro=sum(anthro*fraction), endo=sum(endo*fraction)), by=list(id, continent)]
site_anthro_endo[, human_indoor:= round(anthro*endo, 2)]

setnames(site_anthro_endo, "id", "Site_Name")
site_anthro_endo[, Site_Name:= factor(Site_Name)]
smooth_impact <- merge(smooth_impact, site_anthro_endo, by="Site_Name", all.x=T)
smooth_impact[, human_indoor:=factor(human_indoor)]
write.csv(site_anthro_endo, file.path(ii_dir, "results/tables_for_paper/site_anthro_endo.csv"), row.names = F)

biting_for_plot <- smooth_impact[!Site_Name %in% c(1, 11, 12) & int_id %in% c(101, 25)]
biting_for_plot[, label:= factor(label, labels=c("80% ITN & IRS", "80% AL Case Management"))]

pdf(file.path(out_dir, "biting_intensity.pdf"), width=11, height=7)
ggplot(biting_for_plot, aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.5)+
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=human_indoor, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=human_indoor, group=Site_Name), size=1.25) +
  scale_color_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  scale_fill_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  xlim(0,0.85) +
  ylim(0,0.85) +
  facet_grid(.~label) +
  coord_fixed() +
  theme_minimal(base_size = 8) + 
  labs(x="Initial Prevalence",
       y="Final Prevalence",
       title="")
graphics.off()




### Example Intervention Packages and Maps -----------------------------------------------------------------------------------------------------------------------

ints_to_use <- c(1, 33, 101, 134)

impact_for_plot <- smooth_impact[int_id %in% ints_to_use]

if ("Site_Name" %in% names(model_archs_colormap)){
  impact_for_plot <- merge(impact_for_plot, model_archs_colormap, by="Site_Name", all.x=T)
}else{
  impact_for_plot[, cluster:= as.integer(Site_Name)]
  impact_for_plot <- merge(impact_for_plot, model_archs_colormap, by="cluster", all.x=T)
  
}

impact_for_plot[, cluster_label:= factor(ns_order, labels=model_archs_colormap$name)]
impact_for_plot <- impact_for_plot[cluster %in% sites_to_use]

lineplot_colors <- model_archs_colormap[cluster %in% sites_to_use]$color
lines <- ggplot(impact_for_plot, aes(x=mean_initial, y=mean_final)) +
                geom_abline(size=0.75, alpha=0.25)+
                geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=cluster_label, group=label), alpha=0.25) +
                geom_line(aes(color=cluster_label, linetype=label), size=0.75) +
                scale_color_manual(values=lineplot_colors, name="Site") +
                scale_fill_manual(values=lineplot_colors, name="Site") +
                scale_linetype_manual(values=c( "dotdash", "dotted", "solid","dashed")) + 
                guides(linetype=guide_legend("Intervention"), color = "nonee", fill="none") + 
                xlim(0,0.8) +
                ylim(0,0.8) +
                facet_wrap( . ~ cluster_label ) + 
                theme_minimal(base_size = 8) + 
                theme(legend.position="none",
                      strip.background = element_blank(),
                      # strip.text = element_blank()
                      ) +
                coord_fixed() +
                labs(x="Initial Prevalence",
                     y="Final Prevalence",
                     title="")


impact_dt <- data.table(rasterToPoints(impact_brick[[ints_to_use + 1]]))
impact_dt <- melt(impact_dt, id.vars = c("x", "y"))
setnames(impact_dt, c("x", "y"), c("long", "lat"))
impact_dt[, variable:= factor(variable, labels=c("No Intervention", 
                                                 "40% ITN, 20% IRS, 20% AL CM",
                                                 "80% AL CM",
                                                 "ATSB 3% Kill Rate"))]


color_vals <- generate_full_pal()
prev_cols <- c(color_vals$pal[1:10], rev(brewer.pal(11, "RdYlBu")))
prev_breaks <- c(color_vals$breaks[1:9], seq(0.005, 1, length.out = 12))

maps <- ggplot() +
  geom_raster(data = impact_dt, aes(fill = value, y = lat, x = long)) +
  geom_path(data = africa_dt, aes(x = long, y = lat, group = group), color = "black", size = 0.3) + 
  # scale_fill_gradientn(colors= color_vals$pal, values=color_vals$breaks) +
  scale_fill_gradientn(colors=prev_cols, values=prev_breaks) + 
  facet_wrap(. ~ variable) +  
  coord_equal(xlim = c(-18, 52), ylim = c(-35, 38)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"), 
        legend.title=element_blank(),
        #strip.background = element_blank(),
        #strip.text = element_blank()
        )

lines_height <- ifelse(results_type=="arch", 0.5, 0.35)


pdf(file.path(out_dir, "int_packages.pdf"), width = (8), height = (11))

  lines_vp <- viewport(width = 0.95, height = lines_height, x = 0.45, y = 0.75)
  maps_vp <- viewport(width = 0.95, height = 0.5, x = 0.5, y = 0.25)
  key_vp <- viewport(width = 0.2, height = 0.2, x = 0.85, y = 0.75)
  
  print(lines, vp=lines_vp)
  print(maps, vp = maps_vp)
  print(cluster_plot_for_key, vp=key_vp)
  
graphics.off()


### Megatrends Maps -----------------------------------------------------------------------------------------------------------------------


map_ints_to_use <- c(1, # No interventions
                   105, # 80/80 itn/al_cm
                   125:128, # 80/80/80 + nothing/mAb/pev/tbv
                   132, # 80/80/80 with dp cp
                   134, # 3% atsb + 0/% itn/irs/al_cm
                   152 # 80/80/80 +  25% atsb
        
)

pretty_int_labels <- c("No Intervention", 
                       "80% ITN/AL CM",
                       "80% ITN/AL CM/IRS",
                       "80% ITN/AL CM/IRS + 40% mAB",
                       "80% ITN/AL CM/IRS + 40% TBV",
                       "80% ITN/AL CM/IRS + 40% PEV",
                       "80% ITN/DP CM/IRS",
                       "ATSB 3% Kill Rate",
                       "80% ITN/AL CM/IRS + ATSB 25%"
)

map_impact_dt <- data.table(rasterToPoints(impact_brick[[map_ints_to_use + 1]]))
map_impact_dt <- melt(map_impact_dt, id.vars = c("x", "y"))
setnames(map_impact_dt, c("x", "y"), c("long", "lat"))
map_impact_dt[, variable:= factor(variable, labels=pretty_int_labels)]


color_vals <- generate_full_pal()
prev_cols <- c(color_vals$pal[1:10], rev(brewer.pal(11, "RdYlBu")))
prev_breaks <- c(color_vals$breaks[1:9], seq(0.005, 1, length.out = 12))

maps <- ggplot() +
  geom_raster(data = map_impact_dt, aes(fill = value, y = lat, x = long)) +
  geom_path(data = africa_dt, aes(x = long, y = lat, group = group), color = "black", size = 0.3) + 
  scale_fill_gradientn(colors=prev_cols, values=prev_breaks) + 
  facet_wrap(. ~ variable) +  
  coord_equal(xlim = c(-18, 52), ylim = c(-35, 38)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),            
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"), 
        legend.title=element_blank(),
        # strip.background = element_blank(),
        # strip.text = element_blank()
  )


pdf(file.path(out_dir, "mega_int_packages.pdf"), width = (8), height = (11))
  print(maps)
graphics.off()

### Sensistivity Analysis Figure -----------------------------------------------------------------------------------------------------------------------

if (plot_sens){
  
  
  # plot of all sensitivity sites
  sensitivity_sites <- fread(file.path(sensitivity_dir, "input/site_details.csv" ))
  sensitivity_spoints <- xyFromCell(cluster_raster, sensitivity_sites$id, spatial=T)
  
  pdf(file.path(out_dir, "sensitivity_map.pdf"), width = (5), height = (5))
  
  print(cluster_plot_for_key + 
          geom_point(data=data.table(sensitivity_spoints@coords), aes(x=x, y=y), color="black", shape=3) 
  )
  
  graphics.off()
  
  # load and prep results data
  
  sens_impact <- rbindlist(lapply(c(ii_dir, sensitivity_dir), function(this_dir){
    sites <- fread(file.path(this_dir, "input", "site_details.csv"))
    sites <- sites[continent=="Africa"]
    impact <- fread(file.path(this_dir, "results", "clean", "summary_impact.csv"))
    impact <- smooth_data(impact)
    setnames(impact, "Site_Name", "site_id")
    impact <- merge(sites[, list(site_id=id, cluster)], impact, all.x=T)
    impact[, type:=ifelse(nrow(sites)>50, "sensitivity", "centroid")]
    return(impact)
  }))
  
  # merge colormap
  sens_impact <- merge(sens_impact, model_archs_colormap)
  sens_impact[, cluster_label:= factor(ns_order, labels=model_archs_colormap$name)]
  
  # example plots to show measures of center
  ex_int_to_use <- 79 # 79 is irs: 0%, itn: 60%, al_cm: 60% 
  to_plot <- sens_impact[int_id==ex_int_to_use]
  ex_sens_center_plot <- ggplot(to_plot[type=="sensitivity"], aes(x=mean_initial, y=smooth_mean, group=site_id)) +
    geom_abline() + 
    geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max), alpha=0.25) + 
    geom_line(alpha=0.8) +
    geom_ribbon(data=to_plot[type=="centroid"], aes(ymin=smooth_min, ymax=smooth_max, fill=cluster_label), alpha=0.75, color=NA) + 
    geom_line(data=to_plot[type=="centroid"], aes(color=cluster_label)) +
    scale_color_manual(values=ns_palette) + 
    scale_fill_manual(values=ns_palette) + 
    facet_wrap(~cluster_label) + 
    theme_minimal(base_size = 8) + 
    theme(legend.position = "none") + 
    labs(x="Initial PfPR",
         y="Final PfPR",
         title="")
  
  key_vp <- viewport(width = 0.3, height = 0.3, x = 0.75, y = 0.2)
  
  pdf(file.path(out_dir, paste0("sensitivity_center.pdf")), width=8.5, height=11)
  print (ex_sens_center_plot)
  print(cluster_plot_for_key, vp=key_vp)
  
  graphics.off()
  
  
  # for supplement: full plots to show measures of center
  plot_all_sens <- F
  
  if (plot_all_sens){
    pdf(file.path(out_dir, paste0("supp_sensitivity_center_all.pdf")), width=8.5, height=11)
    
    for (this_id in unique(sens_impact$int_id)){
      print(this_id)
      to_plot <- sens_impact[int_id==this_id]
      
      this_plot <- ggplot(to_plot[type=="sensitivity"], aes(x=mean_initial, y=smooth_mean, group=site_id)) +
        geom_abline() + 
        geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max), alpha=0.25) + 
        geom_line(alpha=0.8) +
        geom_ribbon(data=to_plot[type=="centroid"], aes(ymin=smooth_min, ymax=smooth_max, fill=cluster_label), alpha=0.75, color=NA) + 
        geom_line(data=to_plot[type=="centroid"], aes(color=cluster_label)) +
        scale_color_manual(values=sens_palette) + 
        scale_fill_manual(values=sens_palette) + 
        facet_wrap(~cluster_label) + 
        theme_minimal(base_size = 8) + 
        theme(legend.position = "none") + 
        labs(x="Initial PfPR",
             y="Final PfPR",
             title=unique(to_plot$label))
      print (this_plot)
      print(cluster_plot_for_key, vp=key_vp)
    }
    
    graphics.off()
    
  }
  
  
  
}


###  SVD Vectors -----------------------------------------------------------------------------------------------------------------------

# load svd outputs to get dataset of cell values
load(file.path(arch_dir, "01_svd", "svd_output.rdata")) 

init_variance <- svd_out$d^2/sum(svd_out$d^2)
variance <- data.table(vector=1:length(init_variance), 
                       variance_explained=init_variance)

svd_var_three <- sum(variance[vector<4]$variance_explained)

pdf(file=file.path(out_dir, "svd_variance_explained.pdf"))
varplot <- ggplot(variance[vector<=5], aes(x=vector, y=variance_explained)) +
  geom_line(size=1) +
  geom_point(size=3) +
  theme_minimal(base_size=8) + 
  theme(legend.position = "none") +
  labs(x="Singular Vector", 
       y="Variance Explained",
       title=""
  )
print(varplot)

graphics.off()


###  K-means elbow -----------------------------------------------------------------------------------------------------------------------

this_in_dir <- file.path(arch_dir, "02_kmeans")
cluster_counts <- list.files(this_in_dir)[list.files(this_in_dir) %like% "k_out"]
cluster_counts <- sort(as.integer(gsub("k_out_([0-9]+)_.*", "\\1", cluster_counts)))

# find variance explained
var_explained <- lapply(cluster_counts, function(nclust){
  load(file.path(this_in_dir, paste0("k_out_", nclust, "_cluster", ".rdata")))
  return(k_out$betweenss/k_out$totss) 
})
  

all_var <- data.table(k=cluster_counts, var=unlist(var_explained))
pdf(file=file.path(out_dir, "elbow_plot.pdf"))
elbow_plot<- ggplot(all_var, aes(x=k, y=var)) +
  geom_vline(xintercept = final_nclust, color="blue") + 
  geom_line(size=1) +
  geom_point(shape=1, size=4) +
  theme_minimal(base_size = 8)+
  labs(title= "", 
       x="Cluster Count",
       y="Variance Captured")
print(elbow_plot)
graphics.off()

### Supplement: Covariate Normalization -----------------------------------------------------------------------------------------------------------------------

covariates_final <- melt(svd_wide_datatable, id.vars = c("cov", "variable_name", "variable_val"), variable.name = "id")

covariates_final[, cov:=factor(cov,
                       levels=orig_cov_names,
                       labels=pretty_cov_names)]

cov_type <- ifelse(results_type=="arch", "Rescaled ", "")

pdf(file=file.path(out_dir, "covariate_distributions.pdf"))
  final_covs <- ggplot(covariates_final, aes(x=value)) + 
                geom_density(aes(color=cov, fill=cov), alpha=0.5) +
                facet_grid(cov~., scales="free") +
                theme_minimal() + 
                theme(legend.position = "none") + 
                labs(title= paste0(cov_type, "Distribution of Covariates for SVD"),
                     x= paste0(cov_type,  "Covariate Value"),
                     y="Density") 
  print(final_covs)
graphics.off()

if (results_type=="arch"){
  all_vals <- lapply(orig_cov_names, function(cov_name){
    vals <- fread(file.path(cov_dir, cov_name, paste0(cov_name, "_vals.csv")))
    return(vals)
  })
  
  # keep only those pixels with values for all covariates
  non_null_ids <- lapply(all_vals, function(df){
    return(unique(df$id))
  })
  shared_ids <- Reduce(intersect, non_null_ids)
  all_vals <- rbindlist(all_vals, fill = T)
  all_vals <- all_vals[id %in% shared_ids]
  
  all_vals[variable_name=="month", variable_val:= str_pad(variable_val, 2, side="left", pad="0")]
  all_vals <- all_vals[order(cov, variable_name, variable_val)]
  
  all_vals[, cov:=factor(cov,
                                 levels=orig_cov_names,
                                 labels=pretty_cov_names)]

  
  pdf(file=file.path(out_dir, "unscaled_covariate_distributions.pdf"))
    orig_covs <- ggplot(all_vals, aes(x=value)) + 
      geom_density(aes(color=cov, fill=cov), alpha=0.5) +
      facet_grid(cov~., scales="free") +
      theme_minimal() + 
      theme(legend.position = "none") + 
      labs(title="Unscaled Distribution of Covariates for SVD",
           x="Covariate Value",
           y="Density") 
    print(orig_covs)
  graphics.off()
  
}


### Supplement: All package plots and maps -----------------------------------------------------------------------------------------------------------------------

# takes FOREVER with ggplot framework, only run if necessary

plot_all <- F

if (plot_all){
  
  dir.create(file.path(out_dir, "all_int_plots"), showWarnings = F)
  
  n_perpage <- 4
  start_idx <- 1
  max_int_id <-  max(smooth_impact$int_id)
  
  if ("Site_Name" %in% names(model_archs_colormap)){
    smooth_impact <- merge(smooth_impact, model_archs_colormap, by="Site_Name", all.x=T)
    smooth_impact <- smooth_impact[cluster %in% sites_to_use]
    smooth_impact[, Site_Name:=factor(cluster)]
  }  
  
  while(start_idx < max_int_id){
    end_idx <- start_idx + n_perpage-1
    end_idx <- ifelse( end_idx > max_int_id, max_int_id, end_idx)
    
    print(paste(start_idx, end_idx))
    
    this_impact <- smooth_impact[Site_Name %in% sites_to_use & int_id %in% start_idx:end_idx]
    
    lines <- ggplot(this_impact, aes(x=mean_initial, y=mean_final)) +
      geom_abline(size=1.5, alpha=0.5)+
      geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=Site_Name, group=Site_Name), alpha=0.25) +
      geom_line(aes(color=Site_Name, group=Site_Name), size=1.25) +
      scale_color_manual(values=lineplot_colors, name="Site ID") +
      scale_fill_manual(values=lineplot_colors, name="Site ID") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none") +
      xlim(0,0.85) +
      ylim(0,0.85) +
      facet_grid(factor(label) ~ .) +
      coord_fixed() +
      labs(x="Initial Prevalence",
           y="Final Prevalence",
           title="")
    
    impact_dt <- data.table(rasterToPoints(impact_brick[[(start_idx+1):(end_idx+1)]]))
    names(impact_dt) <- c("long", "lat", as.character(unique(this_impact$label)))
    impact_dt <- melt(impact_dt, id.vars = c("long", "lat"))
    
    maps <- ggplot() +
      geom_raster(data = impact_dt, aes(fill = value, y = lat, x = long)) +
      geom_path(data = africa_dt, aes(x = long, y = lat, group = group), color = "black", size = 0.3) + 
      # scale_fill_gradientn(colors= color_vals$pal, values=color_vals$breaks) +
      scale_fill_gradientn(colors=prev_cols, values=prev_breaks) + 
      facet_grid(variable~ .) +  
      coord_equal(xlim = c(-18, 52), ylim = c(-35, 38)) +
      labs(x = NULL, y = NULL, title = NULL) +
      theme_classic(base_size = 12) +
      theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "in"), legend.title=element_blank())
    
    lines_vp <- viewport(width = 0.4, height = 1, x = 0.4, y = 0.5)
    maps_vp <- viewport(width = 0.4, height = 1, x = 0.8, y = 0.5)
    key_vp <- viewport(width = 0.2, height = 0.2, x = 0.1, y = 0.25)
    
    pdf(file.path(out_dir, paste0("all_int_plots/", "int_", start_idx, "_", end_idx, ".pdf")), width=8.5, height=11)
    
      print(lines, vp=lines_vp)
      print(maps, vp = maps_vp)
      print(cluster_plot_for_key, vp=key_vp)
      
    graphics.off()
    
    start_idx <- start_idx + n_perpage
    
  }
  
  
  
  
  
}




