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

colors <- c("#00a08a", "#d71b5a", "#f2a200", "#902e57", "#f98400","#5392c2")

if (!dir.exists(in_dir)){
  stop("input directory not found! Did you remember to mount your drives?")
}

setwd(in_dir)

raster_to_dt <- function(rast){
  vals <- as.matrix(rast)
  vals <- data.table(id = which(!is.na(vals) & !is.infinite(vals)),
                     value = vals[!is.na(vals) & !is.infinite(vals)])
  vals[, type:= names(rast)]
  return(vals)
}

mask_path <- file.path(out_dir, "mask_layer.tif")
if (!file.exists(mask_path)){
  # set and save a mask for values where access>0 in 2016
  mask_layer <- raster("ITN_2016.ACC.tif")
  mask_layer[mask_layer==0] <- -Inf
  writeRaster(mask_layer, mask_path)
}
mask_layer <- raster(mask_path)


years <- 2010:2016

pdf(file.path(out_dir, "itn_cube_stratification.pdf"), width=10, height=7)
for (year in years){
  print(year)
  basename <- paste0("ITN_", year)
  
  access <- raster::mask(raster(paste0(basename, ".ACC.tif")), mask_layer)
  use <- raster::mask(raster(paste0(basename, ".USE.tif")), mask_layer)
  gap <- raster::mask(raster(paste0(basename, ".GAP.tif")), mask_layer)
  access_gap <- 0.8-access
  use_gap <- access-use
  
  # Sam's description of relative gain
  use_gap_prop <- use_gap/access
  use_gap_prop[access==0] <- 0
  use_gap_prop[use_gap_prop<0] <- 0 # not immediately intuitive what a use gap<0 means
  
  
  # todo: focus on mean and deviation later
  
  stacked_outputs <- stack(access, use, access_gap, use_gap_prop)
  
  nl <- nlayers(stacked_outputs)
  m <- matrix(1:nl, ncol=2, byrow = T)
  layer_names <- c("Access", "Use", "Access Gap (to 80%)", "Use Gap")
  
  for (i in 1:nl){
    index <- which(m==i, arr.ind=T)
    this_rast <- stacked_outputs[[i]]
    
    if (this_rast@data@min<0){
      pal <- wpal("diverging_blue_lightpurple_pink", n=15)[3:13]
      breaks <- c(seq(this_rast@data@min, 0, length.out=8), seq(0, this_rast@data@max, length.out=8)[2:8])
    }else{
      if (layer_names[[i]] %like% "Use Gap"){
        pal <- wpal("diverging_blue_lightpurple_pink", n=15)[8:13]
      }else{
        pal <- wpal("cool_green_grassy", n=15)
      }
      breaks <- seq(0, 1, length.out = 15)
    }
    
    p <- levelplot(this_rast,
                   par.settings=rasterTheme(region=pal), at=breaks,
                   # par.settings=rasterTheme(region=rev(wpal("cool_green_grassy"))),
                   xlab=NULL, ylab=NULL, scales=list(draw=F),
                   main=paste(layer_names[[i]], year), margin=F)
    # print(p, split=c(index[2], index[1], ncol(m), nrow(m)), more=(i<nl))
  }
  
  # play with notion of relative gain
  
  categorical <- copy(access)
  names(categorical) <- "action_class"
  categorical[access>=0.8 & use_gap_prop<=0] <- 1
  categorical[access>=0.8 & use_gap_prop>0] <- 2
  categorical[access<0.8 & use_gap_prop<=0] <- 3
  categorical[access<0.8 & use_gap_prop>0] <- 4
  categorical <- ratify(categorical)
  
  attrs <- levels(categorical)[[1]]
  attrs$action <- c("Do Nothing", "Increase Use", "Increase Access", "Increase Access and Use")
  attrs$category <- c("Access>80%, Use Gap=0 (Do Nothing)",
                      "Access>80%, Use Gap>0 (Increase Use)",
                      "Access<80%, Use Gap=0 (Increase Access)",
                      "Access<80%, Use Gap>0 (Increase Access and Use)")
  levels(categorical) <- attrs
  
  
  catplot <- levelplot(categorical, att="category",
                  col.regions=colors[1:4],
                  xlab=NULL, ylab=NULL, scales=list(draw=F),
                  main=paste("Access/Use Stratification", year), margin=F)
  print(catplot)
  
  
  relgain <- lapply(list(access, use_gap, categorical), raster_to_dt)
  relgain <- rbindlist(relgain)
  relgain <- dcast.data.table(relgain, id ~ type)
  do_both <- relgain[action_class==4]
  names(relgain) <- c("id", "access", "action_class", "use_gap")
  relgain[, action_class:=factor(action_class, labels=attrs$action)]
  
  access_distplot <- ggplot(relgain, aes(x=access, color=action_class, fill=action_class)) +
                    geom_density(alpha=0.75) +
                    facet_grid(action_class~., scales="free") +
                    scale_color_manual(values = colors[1:4]) +
                    scale_fill_manual(values=colors[1:4]) +
                    theme(legend.position = "none") +
                    labs(x="Access",
                         y="Density",
                         title="Distribution of Access by Action Class")
  
  usegap_distplot <- ggplot(relgain[use_gap>0], aes(x=use_gap, color=action_class, fill=action_class)) +
                      geom_density(alpha=0.75) +
                      facet_grid(action_class~., scales="free") +
                      scale_color_manual(values = colors[3:4]) +
                      scale_fill_manual(values=colors[3:4]) +
                      theme(legend.position = "none") +
                      labs(x="Use Gap",
                           y="Density",
                           title="Distribution of Access by Action Class")
  
  
}
graphics.off()







