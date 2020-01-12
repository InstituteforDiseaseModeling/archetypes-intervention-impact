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

## PfPR

metric <- "pfpr"
tifname <- paste0(metric, "_Africa.tif")

new <- brick(file.path(base_dir, new_dir, suffix, tifname))
old <- brick(file.path(base_dir, old_dir, suffix, tifname))
diff <- new-old

nplots <- nlayers(new)

color_vals <- generate_full_pal()

divpal <- c(pnw_palette("Lake", 60)[10:59],  rev(pnw_palette("Shuksan", 35))[3:17], rev(pnw_palette("Starfish", 75))[1:35])
difference_vals <- list(pal=divpal, breaks=seq(-0.5, 0.5, 0.025))
difference_vals_par <- list(pal=divpal, breaks=seq(-0.5, 0.5, 0.025))

# plot_comparison <- function() 

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

pdf(file.path(base_dir, new_dir, suffix, paste0("comparison_", metric, ".pdf")), width=8, height=3)

for (i in 1:nplots){
  new_plot <- levelplot(new[[i]], par.settings=rasterTheme(color_vals$pal), at=color_vals$breaks,
                        xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main="New")
  
  old_plot <- levelplot(old[[i]], par.settings=rasterTheme(color_vals$pal), at=color_vals$breaks,
                        xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main="Old")
  
  diff_plot <- levelplot(diff[[i]], par.settings=rasterTheme(difference_vals$pal), at=difference_vals$breaks,
                         xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main="Diff")
  
  full_plot <- grid.arrange(new_plot, old_plot, diff_plot, ncol=3, top=name_list[i])
  print(full_plot)
  
}


graphics.off()



