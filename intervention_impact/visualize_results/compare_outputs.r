###############################################################################################################
## compare_outputs.r
## Amelia Bertozzi-Villa
## January 2020
## 
##############################################################################################################

library(lattice)
library(data.table)
library(raster)
library(rasterVis)
library(PNWColors)
library(MapSuite)
library(gridExtra)

rm(list=ls())

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results")
setwd(func_dir)
source("pr_to_r0.r")
source("map_ii_functions.r")

old_dir <- "20191008_replicate_megatrends"
new_dir <- "20191009_megatrends_era5_new_archetypes"
suffix <- "results/megatrend_ii_rasters"
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")


## color palettes and breaks
level_palettes <- list(pfpr=generate_full_pal(),
                     par=list(pal=wpal("seaside"),
                              breaks=NA,
                              zscaleLog=T),
                     repro_number=list(pal=c("#e0e0e0", terrain.colors(31)[1:30]),
                                       breaks=c(0,1,seq(1.5, 5, length.out=25), seq(5.1, 80, length.out=4)),
                                       zscaleLog=F))

divpal <- c(pnw_palette("Lake", 60)[10:59],  rev(pnw_palette("Shuksan", 35))[3:17], rev(pnw_palette("Starfish", 75))[1:35])
diff_palettes <- list(pfpr=list(pal=divpal,
                                breaks=seq(-0.5, 0.5, 0.025),
                                zscaleLog=F),
                      par=list(pal=divpal,
                               breaks=(seq(-250000, 250000, 10000)),
                               zscaleLog=F),
                      repro_number=list(pal=divpal,
                                        breaks=seq(-0.5, 0.5, 0.025),
                                        zscaleLog=F)
                      )

# Pretty intervention names
name_list <- c("No Intervention", 
               "80/80 ITN/IRS",
               "80/80/80 ITN/IRS/AL CM",
               "80/80/80 ITN/IRS/DP CM",
               "80/80/80 + 40% mAb",
               "80/80/80 + 40% PEV",
               "80/80/80 + 40% TBV",
               "80/80/80 + 0.15% ATSB",
               "80/80/80 + 3% ATSB",
               "80/80/80 15% ATSB",
               "80/80/80 + 25% % ATSB",
               "3% ATSB + 0/0/0",
               "3% ATSB + 20/20/20",
               "3% ATSB + 40/40/40",
               "3% ATSB +  60/60/60"
)

## make into function

metric <- "par"

for (metric in names(level_palettes)){
  
  print(paste("Comparing for", metric))
  
  tifname <- paste0(metric, "_Africa.tif")
  this_level_pal <- level_palettes[[metric]]
  this_diff_pal <- diff_palettes[[metric]]
  
  new <- brick(file.path(base_dir, new_dir, suffix, tifname))
  old <- brick(file.path(base_dir, old_dir, suffix, tifname))
  diff <- new-old
  
  nplots <- nlayers(new)
  
  print("Plotting")
  pdf(file.path(base_dir, new_dir, suffix, paste0("comparison_", metric, ".pdf")), width=8, height=3)
  
  for (i in 1:nplots){
    
    if (metric=="par"){
      
      new_sum <- cellStats(new[[i]], sum)/1000000
      old_sum <- cellStats(old[[i]], sum)/1000000
      
      new_plot <- levelplot(new[[i]]+1, par.settings=rasterTheme(this_level_pal$pal),
                            zscaleLog=T, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste("New:", round(new_sum, 1), "M"))
      
      old_plot <- levelplot(old[[i]]+1, par.settings=rasterTheme(this_level_pal$pal), 
                            zscaleLog=T, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste("Old:", round(old_sum, 1), "M"))
    }else{
      new_plot <- levelplot(new[[i]], par.settings=rasterTheme(this_level_pal$pal), at=this_level_pal$breaks,
                            zscaleLog=F, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main="New")
      
      old_plot <- levelplot(old[[i]], par.settings=rasterTheme(this_level_pal$pal), at=this_level_pal$breaks,
                            zscaleLog=F, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main="Old")
    }
    
    # diff_plot <- levelplot(abs(diff[[i]])+addition, par.settings=rasterTheme(this_diff_pal$pal), # at=this_diff_pal$breaks,
    #                        zscaleLog=T, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main="Diff")
    
    full_plot <- grid.arrange(new_plot, old_plot, ncol=2, top=name_list[i])
    print(full_plot)
    
  }
  
  graphics.off()
  
 }




