## -----------------------------------------------------------------------------------------------------------------
# ITN Cube
# visualize_cube_results.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# Jan 2019
# 
# Plot and analyze net use, net access, and use gap results from the analysis performed in the "generate_results" 
# scripts in this folder. 
# 
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(rasterVis)
library(gridExtra)
library(MapSuite)

rm(list=ls())

in_dir <- "/Volumes/map_pit/sam/bld1/bras2280/ITNcube"
out_dir <- "/Users/bertozzivill/Dropbox (IDM)/Malaria Team Folder/projects/map_itn_cube"

if (!dir.exists(in_dir)){
  stop("input directory not found! Did you remember to mount your drives?")
}

setwd(in_dir)

mask_path <- file.path(out_dir, "mask_layer.tif")
if (!file.exists(mask_path)){
  # set and save a mask for values where access>0 in 2016
  mask_layer <- raster("ITN_2016.ACC.tif")
  mask_layer[mask_layer==0] <- -Inf
  writeRaster(mask_layer, mask_path)
}
mask_layer <- raster(mask_path)


years <- 2010:2016

pdf(file.path(out_dir, "itn_cube_results.pdf"))
for (year in years){
  print(year)
  basename <- paste0("ITN_", year)
  
  access <- raster::mask(raster(paste0(basename, ".ACC.tif")), mask_layer)
  use <- raster::mask(raster(paste0(basename, ".USE.tif")), mask_layer)
  gap <- raster::mask(raster(paste0(basename, ".GAP.tif")), mask_layer)
  access_gap <- 0.8-access
  # access_gap[access_gap<0] <- -0.001
  testgap <- access-use
  
  # todo: focus on mean and deviation later
  
  stacked_outputs <- stack(access, use, access_gap, testgap)
  
  nl <- nlayers(stacked_outputs)
  m <- matrix(1:nl, ncol=2, byrow = T)
  layer_names <- c("Access", "Use", "Access Gap (to 80%)", "(Access-Use)", "Use Gap?", "Use Gap?")
  
  for (i in 1:nl){
    index <- which(m==i, arr.ind=T)
    this_rast <- stacked_outputs[[i]]
    
    if (this_rast@data@min<0){
      pal <- wpal("diverging_blue_lightpurple_pink", n=15)
      breaks <- c(seq(this_rast@data@min, 0, length.out=8), seq(0, this_rast@data@max, length.out=8)[2:8])
    }else{
      pal <- wpal("cool_green_grassy", n=15)
      breaks <- seq(0, 1, length.out = 15)
    }
    
    p <- levelplot(this_rast,
                   par.settings=rasterTheme(region=pal), at=breaks,
                   # par.settings=rasterTheme(region=rev(wpal("cool_green_grassy"))),
                   xlab=NULL, ylab=NULL, scales=list(draw=F),
                   main=paste(layer_names[[i]], year), margin=F)
    print(p, split=c(index[2], index[1], ncol(m), nrow(m)), more=(i<nl))
  }
  
}
graphics.off()







