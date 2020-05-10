###############################################################################################################
## explore_am_itn_net_count.r
## Amelia Bertozzi-Villa
## April 2020
## 
## See how high other interventions need to be to make up for lowered antimalarial coverage
##############################################################################################################

library(data.table)
library(ggplot2)
library(pracma)
library(raster)
library(rasterVis)

rm(list=ls())

### Setup  ----------------------------------------------------------------------------------------------------------------------------------------------------------

print("Setting up")
func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results")
setwd(func_dir)
source("pr_to_r0.r")
source("map_ii_functions.r")

analysis_subdir <- "20200506_reextract_20191009_mega_era5_new_arch"
analysis_metrics <- c("prev", "severe_inc")
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
archetype_raster_dir <- file.path(base_dir, "../archetypes/results")
main_dir <- file.path(base_dir, analysis_subdir)
suffix <- ""

in_dir <-  file.path(main_dir,"results", "rasters", "base_2020")


all_data <- rbindlist(lapply(analysis_metrics, function(this_metric){
  subset <- fread(file.path(in_dir, this_metric, "itn_cm_scenario_comparisons.csv"))
  return(subset)
}))


all_data[, metric_label :=ifelse(metric=="prev", "Prevalence", "Severe Incidence")]
all_data[, cant_maintain:=factor(cant_maintain, levels=c("Cannot mitigate", "Can mitigate"), labels=c("Cannot mitigate", "Can mitigate"))]

pdf("~/Desktop/compare_all_itn_cm.pdf", width=12, height=8)
for (this_iso in unique(all_data$iso3)){
  print(this_iso)
  this_data <- all_data[iso3==this_iso]
  compare_plot <- ggplot(this_data[scenario %like% "CM"]) +
                  geom_bar(aes(x=scenario, y=nets_in_millions, alpha=scenario, fill=cant_maintain), position="dodge", stat="identity") +
                  geom_hline(data=this_data[!scenario %like% "CM"], aes(yintercept=nets_in_millions, linetype=scenario)) + 
                  facet_grid(.~metric_label) +
                  scale_linetype_manual(values=c("dashed", "solid")) + 
                  scale_alpha_manual(values=c(0.25, 0.5, 0.75, 1)) +
                  theme(axis.text.x = element_text(angle=45, hjust=1),
                        legend.position = "right") +
                  guides(alpha="none", 
                         fill=guide_legend(""),
                         linetype=guide_legend("")) + 
                  labs(x="",
                       y="ITNs (Millions)",
                       title=paste0(this_iso, ": ITNs needed In-Country to Maintain Prevalence and Severe Incidence"))
  
  if (!"Cannot mitigate" %in%  unique(this_data$cant_maintain)){
    compare_plot <- compare_plot + scale_fill_manual(values=c("#00BFC4")) 
  }
  
  
  print(compare_plot)
                
}

graphics.off()






