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
analysis_metric <- "prev"
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
archetype_raster_dir <- file.path(base_dir, "../archetypes/results")
main_dir <- file.path(base_dir, analysis_subdir)
suffix <- ""

in_dir <-  file.path(main_dir,"results", "rasters", "base_2020", analysis_metric)


africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"

itn_main_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20200418_BMGF_ITN_C1.00_R1.00_V2/04_predictions/"
itn_raster_fname <- file.path(itn_main_dir, "ITN_2020_use.tif")
itn_percapita_raster_fname <- file.path(itn_main_dir, "ITN_2020_percapita_nets.tif")
itn_percapita_baseline_raster_fname <- file.path(itn_main_dir, "ITN_2020_percapita_nets.tif")
itn_percapita_no_dists_raster_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20200507_BMGF_ITN_C0.00_R0.00_V2/04_predictions/ITN_2020_percapita_nets.tif"

pop_raster_fname <- "/Volumes/GoogleDrive/Shared drives/MAP Master Outputs/mastergrids/Other_Global_Covariates/Population/Worldpop_GPWv4_Hybrid_201708/5km/Global_Hybrid_Pop_v2_5km_sum_UNAdj_2020.tif"
admin_raster_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/african_cn5km_2013_no_disputes.tif"

# iso-gaul map
iso_gaul_map<-fread("/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/iso_gaul_map.csv")
setnames(iso_gaul_map, c("GAUL_CODE", "COUNTRY_ID", "NAME"), c("gaul", "iso3", "country"))

# by what percent do you want to *reduce* act coverage?
act_percent_reductions <- c(20, 40, 60, 80)

### Shapefile  -----------------------------------------------------
print("Loading Shapefile")
africa_shp <- readOGR(africa_shp_dir)
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)


### Load and format rasters   #####----------------------------------------------------------------------------------------------------------------------------------
print("Loading and formatting rasters")
itn_raster <- raster(itn_raster_fname)

# crop to itn layer, it's the most restricted
itn_raster[itn_raster==0] <- NA


# also find percapita nets, population, and admin units for conversion
itn_percapita_raster <- raster(itn_percapita_raster_fname)
itn_percapita_baseline_raster <- raster(itn_percapita_baseline_raster_fname)
pop_raster <- raster(pop_raster_fname)
pop_raster <- crop(pop_raster, itn_raster)
pop_raster <- raster::mask(pop_raster, itn_raster)

admin_raster <- raster(admin_raster_fname)

# find aggregated net count at baseline
aggregate_npc <- function(npc, admin, admin_map){
  npc_aggregated <- data.table(zonal(npc, admin, fun="sum"))
  names(npc_aggregated) <- c("gaul", "nets")
  npc_aggregated <- merge(npc_aggregated, admin_map, all.x=T)
  npc_aggregated <- npc_aggregated[nets>0]
}

net_count_baseline_raster <- itn_percapita_baseline_raster * pop_raster
net_count_baseline_aggregated <- aggregate_npc(net_count_baseline_raster, admin_raster, iso_gaul_map)
net_count_baseline_aggregated[, type:="Planned Distribution"]

net_count_no_dists_raster <- raster(itn_percapita_no_dists_raster_fname) * pop_raster
net_count_no_dists_aggregated <- aggregate_npc(net_count_no_dists_raster, admin_raster, iso_gaul_map)
net_count_no_dists_aggregated[, type:="No Distribution"]

net_count_baseline_aggregated <- rbind(net_count_baseline_aggregated, net_count_no_dists_aggregated)

### Load counterfactuals, analyze and aggregate  #####----------------------------------------------------------------------------------------------------------------------------------

print("Finding countries that can't maintain prevalence with nets alone")
counterfactuals <- stack(list.files(in_dir, pattern=".tif", full.names=T))

# cap counterfactuals at true baseline coverage
these_names <- names(counterfactuals)
for (idx in 1:nlayers(counterfactuals)){
  counterfactuals[[idx]][itn_raster>counterfactuals[[idx]]] <- itn_raster[itn_raster>counterfactuals[[idx]]]
  # assume that places with true baseline coverage >80% fall into the "can't mitigate" scenario
  counterfactuals[[idx]][counterfactuals[[idx]]>0.8] <- 1
}
names(counterfactuals) <- these_names

cant_maintain <- counterfactuals==1
names(cant_maintain) <- gsub("derived_itn_cov_", "", names(counterfactuals))
names(cant_maintain) <- gsub("act", "cm", names(cant_maintain))

pdf(file.path(in_dir, "cant_mitigate.pdf"), width=7, height=7)
cant_maintain_plot <-levelplot(cant_maintain,
                            par.settings=rasterTheme(region= brewer.pal(9, "Purples")),
                            xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F, main="Pixels where ITN coverage cannot mitigate CM reduction") +
                    latticeExtra::layer(sp.polygons(africa_shp))
print(cant_maintain_plot)
graphics.off()

cant_maintain_pop <- cant_maintain * pop_raster
cant_maintain_pop[[nlayers(cant_maintain_pop)+1]] <- pop_raster
names(cant_maintain_pop) <- c(names(cant_maintain), "pop")

cant_maintain_aggregated <- data.table(zonal(cant_maintain_pop, admin_raster, fun="sum"))
setnames(cant_maintain_aggregated, "zone", "gaul")
cant_maintain_aggregated <- melt(cant_maintain_aggregated, id.vars=c("gaul", "pop"), value.name = "pop_cant_maintain", variable.name = "type")
cant_maintain_aggregated[, type:=paste0("-", gsub(".*_([0-9]+).*", "\\1", type), "% CM")]
cant_maintain_aggregated <- merge(cant_maintain_aggregated, iso_gaul_map, all.x=T)
cant_maintain_aggregated <- cant_maintain_aggregated[pop>0]
setcolorder(cant_maintain_aggregated, c("iso3", "country", "gaul", "type", "pop", "pop_cant_maintain"))
cant_maintain_aggregated <- cant_maintain_aggregated[order(iso3, type)]
cant_maintain_aggregated[, perc_cant_maintain:=pop_cant_maintain/pop]

maintain_cutoff <- 0.25
cant_maintain_aggregated[, cant_maintain_binary:=ifelse(perc_cant_maintain>=maintain_cutoff, "Cannot maintain", "Can maintain")]


print("Calculating nets needed")
counterfactuals_aggregated <- lapply(1:nlayers(counterfactuals), function(raster_idx){
  
  this_itn_raster <- counterfactuals[[raster_idx]]
  percent_reduction <- as.integer(gsub(".*_([0-9]+).*", "\\1", names(this_itn_raster)))
  print(percent_reduction)
  
  perc_increase_itn <- this_itn_raster-itn_raster
  perc_increase_itn[perc_increase_itn<0] <- 0
  
  perc_change_itn <- this_itn_raster/itn_raster
  perc_change_itn[perc_change_itn<1] <- 1
  this_percapita_itn <- itn_percapita_raster*perc_change_itn
  this_percapita_itn[this_percapita_itn>1] <- 1
  
  
  this_net_count <- this_percapita_itn * pop_raster
  this_net_count_aggregated <- aggregate_npc(this_net_count, admin_raster, iso_gaul_map)
  this_net_count_aggregated[, type:=paste0("-", percent_reduction, "% CM")]
  
  return(this_net_count_aggregated)
})

counterfactuals_aggregated <- rbindlist(counterfactuals_aggregated)
compare_aggregated <- rbind(net_count_baseline_aggregated, counterfactuals_aggregated)
compare_aggregated[, class:=ifelse(type %like% "CM", "counter", "baseline")]
compare_aggregated[, mil_nets:=nets/1000000]

compare_aggregated <- merge(compare_aggregated, cant_maintain_aggregated[, list(iso3, type, perc_cant_maintain, cant_maintain=cant_maintain_binary)], all.x=T)
# compare_aggregated[class=="baseline", cant_maintain:=0]
compare_aggregated[, cant_maintain:=factor(cant_maintain, levels=c("Cannot maintain", "Can maintain"), labels=c("Cannot mitigate", "Can mitigate"))]

compare_aggregated[iso3=="CAF", country:="CAR"]
compare_aggregated[iso3=="COD", country:="DRC"]

label <- ifelse(analysis_metric=="prev", "Prevalence", "Severe Incidence")

pdf(file.path(in_dir, "mitigation_bar_plots.pdf"), width=14, height=10)
ggplot(compare_aggregated[class=="counter"]) +
  geom_bar(aes(x=type, y=mil_nets, alpha=type, fill=cant_maintain), position="dodge", stat="identity") +
  geom_hline(data=compare_aggregated[class=="baseline"], aes(yintercept=mil_nets, linetype=type)) + 
  facet_wrap(~country, scales="free_y") +
  scale_linetype_manual(values=c("dashed", "solid")) + 
  # scale_fill_manual(values=c("#00BFC4")) + 
  scale_alpha_manual(values=c(0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "right") +
  guides(alpha="none", 
         fill=guide_legend(""),
         linetype=guide_legend("")) + 
  labs(x="",
       y="ITNs (Millions)",
       title=paste("ITNs needed In-Country to Mitigate", label))
graphics.off()

# compare <- stack(itn_percapita_raster, this_percapita_itn)
# names(compare) <- c("NPC 2019", "NPC Needed")
# 

for_bmgf <- compare_aggregated[, list(iso3, country, metric=analysis_metric, scenario=type, 
                                      nets, nets_in_millions=mil_nets,
                                      perc_cant_maintain, cant_maintain)]
write.csv(for_bmgf, file.path(in_dir, "itn_cm_scenario_comparisons.csv"), row.names=F)



