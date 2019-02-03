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
zed_root <- '/Volumes/map_data/'
compare_dir <- file.path(zed_root, 'cubes/5km/ITN/')
out_dir <- "/Users/bertozzivill/Dropbox (IDM)/Malaria Team Folder/projects/map_itn_cube"

colors <- c("#00a08a", "#d71b5a", "#f2a200", "#902e57", "#f98400","#5392c2")

if (!dir.exists(in_dir)){
  stop("input directory not found! Did you remember to mount your drives?")
}


years <- 2000:2016

pdf(file.path(out_dir, "use_compare.pdf"), width=7, height=10)
for (year in years){
  print(year)
  basename <- paste0("ITN_", year)
  sam_use <- raster(file.path(in_dir, paste0(basename, ".USE.tif")))
  zed_use <- raster(file.path(compare_dir, paste0(year, '.ITN.use.yearavg.new.adj.tif')))
  diff <- sam_use-zed_use
  
  stacked_outputs <- stack(sam_use, zed_use)
  names(stacked_outputs) <- c("Sam Dir", "Z Dir")
  
  pal <- brewer.pal(8, "RdYlGn")
  breaks <- seq(0, 1, length.out = length(pal)+1)
  
  p <- levelplot(stacked_outputs,
            par.settings=rasterTheme(region=pal), at=breaks,
            xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F,
            main=paste(year, ": Min Diff", diff@data@min, ", Max Diff", diff@data@max)
            )
  print(p)
  
}
graphics.off()

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

pdf(file.path(out_dir, "access_vs_use.pdf"), width=7, height=10)
for (year in years){
  print(year)
  basename <- paste0("ITN_", year)
  
  access <- raster::mask(raster(paste0(basename, ".ACC.tif")), mask_layer)
  use <- raster::mask(raster(paste0(basename, ".USE.tif")), mask_layer)
  gap <- raster::mask(raster(paste0(basename, ".GAP.tif")), mask_layer)
  access_gap <- 0.8-access
  use_gap <- access-use
  
  # use gap as a proportion
  use_gap_prop <- use_gap/access
  use_gap_prop[access==0] <- 0
  use_gap_prop[use_gap_prop<0] <- 0 # not immediately intuitive what a use gap<0 means
  
  # my definitions of access gain and use gain:
  
  use[use>0 & use<1.0e-5] <- 1.0e-5
  capped_use_perc <- min(use, access)*100
  access_perc <- access*100
  access_thresh <- 80
  use_rate <- capped_use_perc/access_perc
  access_diff <- use_rate * access_thresh
  access_diff[access_perc >= access_thresh] <- use_rate[access_perc >= access_thresh]*100
  access_gain <- access_diff-capped_use_perc
  use_gain <- access_perc-capped_use_perc
  
  categorical_gain <- access_gain>use_gain
  categorical_gain <- ratify(categorical_gain)
  attrs <- levels(categorical_gain)[[1]]
  attrs$action <- c("Increase Use", "Increase Access")
  levels(categorical_gain) <- attrs
  
  stacked_outputs <- stack(capped_use_perc, categorical_gain,
                           access_diff, access_gain,
                           access_perc, use_gain)
  
  nl <- nlayers(stacked_outputs)
  m <- matrix(1:nl, ncol=2, byrow = T)
  layer_names <- c("True Use", "Action Needed",
                   "Maximize Access (to 80%)", "Increase in Access Needed",
                   "Maximize Use", "Increase in Use Needed")
  
  for (i in 1:nl){
    index <- which(m==i, arr.ind=T)
    this_rast <- stacked_outputs[[i]]
    this_name <- layer_names[[i]]
    this_name <- ifelse(this_name=="True Use", paste(this_name, year), this_name)
    
    if (this_name %like% "Action"){
      p <- levelplot(this_rast, att="action",
                col.regions=colors[c(1,3)],
                xlab=NULL, ylab=NULL, scales=list(draw=F),
                main=this_name, margin=F)
    }else{
      
      if (this_name %like% "True Use" | this_name %like% "Maximize"){
        pal <- wpal("cool_green_grassy", n=15)
      }else{
        pal <- brewer.pal(9, "RdPu")
      }
      
      p <- levelplot(this_rast,
                     par.settings=rasterTheme(region=pal), at=seq(0, 100, length.out = length(pal)+1),
                     xlab=NULL, ylab=NULL, scales=list(draw=F),
                     main=paste(this_name), margin=F)
      
      
    }
    
    
    print(p, split=c(index[2], index[1], ncol(m), nrow(m)), more=(i<nl))
    
  }

  

}
graphics.off()

# quick compare: national access means and access deviation

mean_access_00 <- raster::mask(raster(paste0("ITN_2000", ".MEAN.tif")), mask_layer)
access_dev_00 <- raster::mask(raster(paste0("ITN_2000", ".DEV.tif")), mask_layer)

mean_access_16 <- raster::mask(raster(paste0("ITN_2016", ".MEAN.tif")), mask_layer)
access_dev_16 <- raster::mask(raster(paste0("ITN_2016", ".DEV.tif")), mask_layer)

stacked_outputs <- stack(mean_access_00, access_dev_00,
                         mean_access_16, access_dev_16)

years <- 2000:2016
nl <- length(years)
m <- matrix(1:20, nrow=4, byrow = T)

for (i in 1:nl){
  
  year <- years[[i]]
  index <- which(m==i, arr.ind=T)
  this_rast <- raster::mask(raster(paste0("ITN_", year, ".MEAN.tif")), mask_layer)
  this_name <- year

  if (this_name %like% "Mean Access" | this_name==year){
    pal <- brewer.pal(8, "RdYlGn")
    breaks <- seq(0, mean_access_00@data@max, length.out = length(pal)+1)
  }else{
    pal <- brewer.pal(9, "PRGn")
    seg_length <- ceiling(length(pal)/2)
    breaks <-  c(seq(-4, 0, length.out = seg_length),
                 seq(0, 4, length.out = seg_length)[2:seg_length])
  }
  
  p <- levelplot(this_rast,
                 par.settings=rasterTheme(region=pal), at=breaks,
                 xlab=NULL, ylab=NULL, scales=list(draw=F),
                 main=paste(this_name), margin=F)
  
  
  print(p, split=c(index[2], index[1], ncol(m), nrow(m)), more=(i<nl))

}

