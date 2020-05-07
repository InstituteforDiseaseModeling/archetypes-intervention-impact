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
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
archetype_raster_dir <- file.path(base_dir, "../archetypes/results")
main_dir <- file.path(base_dir, analysis_subdir)
suffix <- ""
in_dir <-  file.path(main_dir,"results", "rasters")


africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"

itn_main_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20200418_BMGF_ITN_C1.00_R1.00_V2/04_predictions/"
itn_raster_fname <- file.path(itn_main_dir, "ITN_2019_use.tif")
itn_percapita_raster_fname <- file.path(itn_main_dir, "ITN_2019_percapita_nets.tif")
itn_percapita_baseline_raster_fname <- file.path(itn_main_dir, "ITN_2020_percapita_nets.tif")

pop_raster_fname <- "/Volumes/GoogleDrive/Shared drives/MAP Master Outputs/mastergrids/Other_Global_Covariates/Population/Worldpop_GPWv4_Hybrid_201708/5km/Global_Hybrid_Pop_v2_5km_sum_UNAdj_2020.tif"
admin_raster_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/african_cn5km_2013_no_disputes.tif"

# iso-gaul map
iso_gaul_map<-fread("/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/iso_gaul_map.csv")
setnames(iso_gaul_map, c("GAUL_CODE", "COUNTRY_ID", "NAME"), c("gaul", "iso3", "country"))

# by what percent do you want to *reduce* act coverage?
act_percent_reductions <- c(20, 40, 60, 80)

### Load and format rasters  #####----------------------------------------------------------------------------------------------------------------------------------
print("Loading and Formatting rasters")
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
net_count_baseline_raster <- itn_percapita_baseline_raster * pop_raster
net_count_baseline_aggregated <- data.table(zonal(net_count_baseline_raster, admin_raster, fun="sum"))
names(net_count_baseline_aggregated) <- c("gaul", "nets")
net_count_baseline_aggregated <- merge(net_count_baseline_aggregated, iso_gaul_map, all.x=T)
net_count_baseline_aggregated <- net_count_baseline_aggregated[nets>0]
net_count_baseline_aggregated[, type:="Baseline 2020"]


### Load counterfactuals, analyze and aggregate  #####----------------------------------------------------------------------------------------------------------------------------------

print("Finding countries that can't maintain prevalence with nets alone")
counterfactuals <- stack(list.files(in_dir, pattern=".tif", full.names=T))

cant_maintain <- counterfactuals==1
names(cant_maintain) <- names(counterfactuals)

cant_maintain_pop <- cant_maintain * pop_raster
cant_maintain_pop[[nlayers(cant_maintain_pop)+1]] <- pop_raster
names(cant_maintain_pop) <- c(names(counterfactuals), "pop")

cant_maintain_aggregated <- data.table(zonal(cant_maintain_pop, admin_raster, fun="sum"))
setnames(cant_maintain_aggregated, "zone", "gaul")
cant_maintain_aggregated <- melt(cant_maintain_aggregated, id.vars=c("gaul", "pop"), value.name = "pop_cant_maintain", variable.name = "type")
cant_maintain_aggregated[, type:=paste0("-", gsub(".*_([0-9]+)", "\\1", type), "% ACT")]
cant_maintain_aggregated <- merge(cant_maintain_aggregated, iso_gaul_map, all.x=T)
cant_maintain_aggregated <- cant_maintain_aggregated[pop>0]
setcolorder(cant_maintain_aggregated, c("iso3", "country", "gaul", "type", "pop", "pop_cant_maintain"))
cant_maintain_aggregated <- cant_maintain_aggregated[order(iso3, type)]
cant_maintain_aggregated[, perc_cant_maintain:=pop_cant_maintain/pop]

maintain_cutoff <- 0.5
cant_maintain_aggregated[, cant_maintain_binary:=ifelse(perc_cant_maintain>=maintain_cutoff, 1, 0)]


print("Calculating nets needed")
counterfactuals_aggregated <- lapply(1:nlayers(counterfactuals), function(raster_idx){
  
  this_itn_raster <- counterfactuals[[raster_idx]]
  percent_reduction <- as.integer(gsub(".*_([0-9]+)", "\\1", names(this_itn_raster)))
  print(percent_reduction)
  
  perc_increase_itn <- this_itn_raster-itn_raster
  perc_increase_itn[perc_increase_itn<0] <- 0
  
  perc_change_itn <- this_itn_raster/itn_raster
  this_percapita_itn <- itn_percapita_raster*perc_change_itn
  this_percapita_itn[this_percapita_itn>1] <- 1
  
  
  this_net_count <- this_percapita_itn * pop_raster
  this_net_count_aggregated <- data.table(zonal(this_net_count, admin_raster, fun="sum"))
  names(this_net_count_aggregated) <- c("gaul", "nets")
  this_net_count_aggregated <- merge(this_net_count_aggregated, iso_gaul_map, all.x=T)
  this_net_count_aggregated <- this_net_count_aggregated[nets>0]
  this_net_count_aggregated[, type:=paste0(percent_reduction, "% ACT reduction")]
  
  return(this_net_count_aggregated)
})

counterfactuals_aggregated <- rbindlist(counterfactuals_aggregated)
compare_aggregated <- rbind(net_count_baseline_aggregated, counterfactuals_aggregated)
compare_aggregated[, type:=gsub(" reduction", "", type)]
compare_aggregated[type %like% "ACT", type:=paste0("-", type)]
compare_aggregated[type %like% "Baseline", type:="Planned"]
compare_aggregated[, class:=ifelse(type %like% "Planned", "baseline", "counter")]
compare_aggregated[, mil_nets:=nets/1000000]
factor_names <- c("Planned", 
                  "-20% ACT",
                  "-40% ACT",
                  "-60% ACT",
                  "-80% ACT")
compare_aggregated[, type:=factor(type, levels=factor_names, labels=factor_names)]

compare_aggregated <- merge(compare_aggregated, cant_maintain_aggregated[, list(iso3, type, cant_maintain=cant_maintain_binary)], all.x=T)
compare_aggregated[class=="baseline", cant_maintain:=0]
compare_aggregated[, cant_maintain:=factor(cant_maintain, levels=c(1,0), labels=c(1,0))]

ggplot(compare_aggregated[class=="counter"], aes(x=type, y=mil_nets)) +
  geom_bar(aes(alpha=type, fill=cant_maintain), position="dodge", stat="identity") +
  geom_hline(data=compare_aggregated[class=="baseline"], aes(yintercept=mil_nets)) + 
  facet_wrap(~country, scales="free_y") +
  # scale_fill_manual(values=c("#B6ECEE",
  #                            "#91E3E5",
  #                            "#48D1D4", 
  #                            "#00BFC4")) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none") +
  labs(x="",
       y="LLINs to Maintain Prevalence (Millions)")

compare <- stack(itn_percapita_raster, this_percapita_itn)
names(compare) <- c("NPC 2019", "NPC Needed")


