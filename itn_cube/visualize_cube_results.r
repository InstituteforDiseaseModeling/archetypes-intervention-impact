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
  access_diff <- (capped_use_perc/access_perc) * 100
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
                   "Maximize Access", "Increase in Access Needed",
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

  
  

  # # todo: focus on mean and deviation later
  # 
  # # play with notion of relative gain
  # 
  # categorical <- copy(access)
  # names(categorical) <- "action_class"
  # categorical[access>=0.8 & use_gap_prop<=0] <- 1
  # categorical[access>=0.8 & use_gap_prop>0] <- 2
  # categorical[access<0.8 & use_gap_prop<=0] <- 3
  # categorical[access<0.8 & use_gap_prop>0] <- 4
  # categorical <- ratify(categorical)
  # 
  # attrs <- levels(categorical)[[1]]
  # attrs$action <- c("Do Nothing", "Increase Use", "Increase Access", "Increase Access and Use")
  # attrs$category <- c("Access>80%, Use Gap=0 (Do Nothing)",
  #                     "Access>80%, Use Gap>0 (Increase Use)",
  #                     "Access<80%, Use Gap=0 (Increase Access)",
  #                     "Access<80%, Use Gap>0 (Increase Access and Use)")
  # levels(categorical) <- attrs
  # 
  # 
  # catplot <- levelplot(categorical, att="category",
  #                 col.regions=colors[1:4],
  #                 xlab=NULL, ylab=NULL, scales=list(draw=F),
  #                 main=paste("Access/Use Stratification", year), margin=F)
  # # print(catplot)
  # 
  # 
  # relgain <- lapply(list(access, use_gap, categorical), raster_to_dt)
  # relgain <- rbindlist(relgain)
  # relgain <- dcast.data.table(relgain, id ~ type)
  # do_both <- relgain[action_class==4]
  # names(relgain) <- c("id", "access", "action_class", "use_gap")
  # relgain[, action_class:=factor(action_class, labels=attrs$action)]
  # 
  # access_distplot <- ggplot(relgain, aes(x=access, color=action_class, fill=action_class)) +
  #                   geom_density(alpha=0.75) +
  #                   facet_grid(action_class~., scales="free") +
  #                   scale_color_manual(values = colors[1:4]) +
  #                   scale_fill_manual(values=colors[1:4]) +
  #                   theme(legend.position = "none") +
  #                   labs(x="Access",
  #                        y="Density",
  #                        title="Distribution of Access by Action Class")
  # 
  # usegap_distplot <- ggplot(relgain[use_gap>0], aes(x=use_gap, color=action_class, fill=action_class)) +
  #                     geom_density(alpha=0.75) +
  #                     facet_grid(action_class~., scales="free") +
  #                     scale_color_manual(values = colors[3:4]) +
  #                     scale_fill_manual(values=colors[3:4]) +
  #                     theme(legend.position = "none") +
  #                     labs(x="Use Gap",
  #                          y="Density",
  #                          title="Distribution of Access by Action Class")
  # 
  # 
}
graphics.off()







