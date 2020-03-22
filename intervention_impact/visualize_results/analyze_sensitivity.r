###############################################################################################################
## analyze_sensitivity.r
## Amelia Bertozzi-Villa
## March 2020
## 
## 
##############################################################################################################


library(data.table)
library(ggplot2)

rm(list=ls())

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")

centroid_dir <- file.path(main_dir, "20191009_mega_era5_new_arch")
sensitivity_dir <- file.path(main_dir, "20191218_site_sensitivity")

full_impact <- rbindlist(lapply(c(centroid_dir, sensitivity_dir), function(this_dir){
  sites <- fread(file.path(this_dir, "input", "site_details.csv"))
  sites <- sites[continent=="Africa"]
  impact <- fread(file.path(this_dir, "results", "clean", "summary_impact.csv"))
  setnames(impact, "Site_Name", "site_id")
  impact <- merge(sites[, list(site_id=id, cluster)], impact, all.x=T)
  impact[, type:=ifelse(nrow(sites)>50, "sensitivity", "centroid")]
  return(impact)
}))


colormap <- fread(file.path(sensitivity_dir, "results", "clean", "cluster_color_map.csv"))
colormap <- colormap[order(ns_order)]
palette <- colormap$color
full_impact <- merge(full_impact, colormap)
full_impact[, cluster_label:= factor(ns_order, labels=colormap$name)]

pdf(file.path(sensitivity_dir, "results", "clean", "sensitivity_results.pdf"), width=8, height=8)
  
  for (this_id in unique(full_impact$int_id)){
    print(this_id)
    to_plot <- full_impact[int_id==this_id]
    this_plot <- ggplot(to_plot[type=="sensitivity"], aes(x=mean_initial, y=mean_final, group=site_id)) +
      geom_abline() + 
      geom_ribbon(aes(ymin=min_final, ymax=max_final), alpha=0.25) + 
      geom_line(alpha=0.8) +
      geom_ribbon(data=to_plot[type=="centroid"], aes(ymin=min_final, ymax=max_final, fill=cluster_label), alpha=0.75, color=NA) + 
      geom_line(data=to_plot[type=="centroid"], aes(color=cluster_label)) +
      scale_color_manual(values=palette) + 
      scale_fill_manual(values=palette) + 
      facet_wrap(~cluster_label) + 
      theme_minimal() + 
      theme(legend.position = "none") + 
      labs(x="Initial PfPR",
           y="Final PfPR",
           title=unique(to_plot$label))
    print (this_plot)
  }

graphics.off()



pdf(file.path(sensitivity_dir, "results", "clean", "sensitivity_distribution.pdf"), width=10, height=8)

for (this_id in unique(full_impact$int_id)){
  print(this_id)
  to_plot <- full_impact[int_id==this_id]
  
  to_plot_sensitivity <- to_plot[type=="sensitivity"]
  setnames(to_plot_sensitivity, "cluster_label", "sensitivity_label")
  to_plot_centroid <- rbindlist(lapply(unique(to_plot$cluster_label), function(this_label){
    subset <- to_plot[type=="centroid"]
    subset[, sensitivity_label:=this_label]
    return(subset)
  }))
  
  this_plot <- ggplot(to_plot_centroid, aes(x=mean_initial, y=mean_final, group=site_id)) +
    geom_abline() + 
    geom_ribbon(aes(ymin=min_final, ymax=max_final, fill=cluster_label), alpha=0.3, color=NA) + 
    geom_line(aes(color=cluster_label)) + 
    geom_ribbon(data=to_plot_sensitivity, aes(ymin=min_final, ymax=max_final), alpha=0.75) + 
    geom_line(data=to_plot_sensitivity) + 
    scale_color_manual(values=palette) + 
    scale_fill_manual(values=palette) + 
    facet_wrap(~sensitivity_label) + 
    theme_minimal() + 
    theme(legend.title = element_blank()) + 
    labs(x="Initial PfPR",
         y="Final PfPR",
         title=unique(to_plot$label))
  print(this_plot)
}

graphics.off()









