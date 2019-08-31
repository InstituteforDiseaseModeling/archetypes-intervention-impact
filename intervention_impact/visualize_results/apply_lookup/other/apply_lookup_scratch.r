
library(data.table)
library(raster)
library(rasterVis)
library(dichromat)
library(MapSuite) # cubehelix palettes; for installation instructions see https://rdrr.io/github/RebeccaStubbs/MapSuite

rm(list=ls())

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results/apply_lookup")
setwd(func_dir)
source("pr_to_r0.r")

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))

main_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/")
lookup_dir <- file.path(main_dir, 
                        "lookup_tables/interactions")
orig_pr_dir <- file.path(main_dir, 
                         "megatrends/pete_analysis/actual_ssp2_base2016_2050.tif")
africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/shapefiles/Africa.shp"
# mask_dir <- "/Volumes/GoogleDrive/Shared drives/MAP Master Outputs/master_geometries/Cartographic/Land_Sea_Masks/admin2018_sea_mask.shp"

# out_dir <- file.path(main_dir, "writing_and_presentations/megatrends/megatrends_gates_2019")
base_label <- "Megatrends Base 2016"
interventions <- c("Baseline:0%, ATSB Initial Kill:1%" # ,
                   # "Baseline:20%, ATSB Initial Kill:1%",
                   # "Baseline:40%, ATSB Initial Kill:1%",
                   # "Baseline:60%, ATSB Initial Kill:1%"
                   ) 

# set 'repro' to true if you want a reproductive number rather than a pfpr estimate
repro <- F
if (repro){
  print("getting repro spline")
  repro_spline <- R2spline()
}

africa_shp <- readOGR(africa_shp_dir)
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)


#################################################################################################
get_splinefunction<-function(X,Y){
  
  plot(X,Y,type="b")
  
  splObj<-smooth.spline(X, Y)
  
  Xpred<-seq(0,1,length=100)
  Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
  # lines(Xpred,Ypred,col=2)
  return(splObj)
}
##################################################################################################

##################################################################################################
apply_spline_to_raster<-function(splObj,inRaster, return_r0=F){
  
  inMat<-as.matrix(inRaster)
  inVec<-as.vector(inMat)
  NAid<-is.na(inVec)
  inVecNoNA<-inVec[!NAid]
  
  outVecNoNA<-predict(splObj,inVecNoNA)$y
  if (return_r0==F){
    outVecNoNA <- pmin(outVecNoNA, inVecNoNA)
    outVecNoNA[outVecNoNA<0] <- 0
  }
  
  outVec<-rep(NA,length(inVec))
  outVec[!NAid]<-outVecNoNA
  outMat<-matrix(outVec,nrow=nrow(inMat))
  outRaster<-inRaster
  outRaster[]<-outMat
  
  return(outRaster)
}
##################################################################################################

# include atsb runs 
atsb_runs <- c("MAP_For_Symposium_ATSB_Higher_Existing_Intervention.csv", 
               "MAP_For_Symposium_ATSB_Lower_Intervention.csv",
               "MAP_For_Symposium_ATSB_Lower_Existing_Intervention.csv",
               "MAP_For_Symposium_ATSB_No_Existing_Intervention.csv")

initial <- fread(file.path(lookup_dir, "../initial/MAP_II_New_Sites_Burnin.csv"))
prelim_data <- rbindlist(lapply(atsb_runs, function(fname){fread(file.path(lookup_dir, fname))}))

lut <- merge(prelim_data, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
lut[, Run_Number:=factor(Run_Number)]
lut[, Intervention:= paste0("Baseline:", ITN_Coverage*100, "%, ", "ATSB Initial Kill:", ATSB_Initial_Effect*100, "%")]

lut[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]
lut[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]

lut <- lut[Intervention %in% interventions]

# read in rasters: PR, cluster assignments, and continent masks
pr_orig <- raster(orig_pr_dir)

cluster_map <- raster("africa_clusters_v4.tif")
# africa counterfactual is smaller than cluster map; align extents and mask oceans etc.
cluster_map <- crop(cluster_map, pr_orig)
pr_orig <- crop(pr_orig, cluster_map)
pr_orig <- raster::mask(pr_orig, cluster_map)


bounding_pr <- raster("actual_ssp2_base2016_2050.tif")
bounding_pr <- crop(bounding_pr, cluster_map)
bounding_pr <- raster::mask(bounding_pr, cluster_map)

continent <- "africa"
cluster_list <-  list("aba"=1,
                      "kasama"=3,
                      "djibo"=4,
                      "kananga"=2,
                      "moine"=6,
                      "gode"=5)


# I want the first panel of my plot to be intervention-delected values of either pfpr or r0
if (repro==T){
  print("finding r0 for raster")
  r0 <-apply_spline_to_raster(repro_spline, pr_orig, return_r0=T)
  raster_list <- list(r0)
}else{
  raster_list <- list(pr_orig)
}

idx <- 2

for (intervention in interventions){
  print(intervention)
  
  # loop through each cluster, apply the appropriate spline, convert to reproductive number if desired
  final_prs <- lapply(names(cluster_list), function(site_name){
    this_mask <- cluster_map==cluster_list[[site_name]]
    print(paste(site_name, ":", cluster_list[[site_name]]))
    pr_masked <- raster::mask(pr_orig, this_mask, maskvalue=0)
    this_lut <- unique(lut[Site_Name==site_name & Intervention==intervention, list(mean_initial, mean_final)])
    this_spline <- get_splinefunction(this_lut$mean_initial,this_lut$mean_final)
    this_pr <- apply_spline_to_raster(this_spline, pr_masked, return_r0 = F)
    
    if (repro==T){
      this_repro <- apply_spline_to_raster(repro_spline, this_pr, return_r0 = T)
      return(this_repro)
    }else{
      return(this_pr)
    }
    
  })
  
  # merge back into a single raster
  this_pr_final <- do.call(merge, final_prs)
  this_pr_final <- min(stack(this_pr_final, bounding_pr)) 
  # this_pr_final <- mask(this_pr_final, cluster_map, maskvalue=FALSE)
  raster_list[[idx]] <- this_pr_final
  idx <- idx + 1 
}

stacked_layers <- stack(raster_list)

stacked_layers <- extend(stacked_layers, africa_shp)
# save and plot

map_pal <-c("#E6E6E6", "#A3A3A3", "#5A17FD", "#51C8FE", "#AEFEAB", "#FBD817",  "#FE3919") 
map_breaks <- c(0, 0.005, 0.01, 0.25, 0.5, 0.75, 1)

generate_ramp <- function(colors, vals, n.out=10){
  rampfunc <- colorRampPalette(colors)
  pal <- rampfunc(n.out)
  breaks <- seq(vals[1], vals[2], length.out=n.out)
  return(list(pal, breaks))
}

full_pal <- NULL
full_breaks <- NULL
n.out=10

for (idx in 1: (length(map_pal)-1) ){
  colors <- map_pal[idx:(idx+1)]
  vals <- map_breaks[idx:(idx+1)]
  ramps <- generate_ramp(colors, vals, n.out)
  
  if (idx == (length(map_pal)-1)){
    full_pal[[idx]] <- ramps[[1]]
    full_breaks[[idx]] <- ramps[[2]]
  }else{
    full_pal[[idx]] <- ramps[[1]][1:(n.out-1)]
    full_breaks[[idx]] <- ramps[[2]][1:(n.out-1)]
  }
} 
full_pal <- unlist(full_pal)
full_breaks <- unlist(full_breaks)

if (repro==T){
  names(stacked_layers) <- c("R0 2015", interventions)
  breaks <- c(seq(-0.001, 5, length.out=50), seq(5, 10, length.out = 46)[2:46], seq(10, 75, length.out=6)[2:6])
  pdf(file.path(out_dir, paste0("repro_", continent, ".pdf")), width=12, height=6)
  print(levelplot(stacked_layers, par.settings=rasterTheme(region=wpal("sky")), zscaleLog=T, xlab=NULL, ylab=NULL, scales=list(draw=F)))
  graphics.off()
  writeRaster(stacked_layers, options="INTERLEAVE=BAND", file.path(out_dir, paste0("repro_",continent, ".tif")), overwrite=T)
}else{
  names(stacked_layers) <- c(base_label, interventions)
  pal <- c("#e0e0e0", rev(brewer.pal(11, "Spectral")))
  breaks <- c(0, seq(0.01, 1, length.out=11))
  pdf(file.path(out_dir, paste0("atsb_pfpr_", continent, ".pdf")), width=12, height=6)
  print(levelplot(stacked_layers[[1]], par.settings=rasterTheme(full_pal), at=full_breaks, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F)) +
          latticeExtra::layer(sp.polygons(africa_shp)) )
  print(levelplot(stacked_layers[[2:5]], par.settings=rasterTheme(full_pal), at=full_breaks, xlab=NULL, margin=F, ylab=NULL, scales=list(draw=F)) +
          latticeExtra::layer(sp.polygons(africa_shp)))
  graphics.off()
  writeRaster(stacked_layers, options="INTERLEAVE=BAND", file.path(out_dir, paste0("abv_pfpr_",continent, ".tif")), overwrite=T, bylayer=T, suffix="names")
}

# plot population at risk and case counts
pal <- wpal("seaside")
pop <- raster("other/ssp5_total_2050_MG_5K.tif")
pop <- crop(pop, stacked_layers)

# has_transmission <- copy(stacked_layers)
# has_transmission[has_transmission<0.01] <- 0
# has_transmission[has_transmission>0] <- 1
# don't know why arithmetic on rasterstacks doesn't work right now, but:
has_transmission <- lapply(1:nlayers(stacked_layers), function(idx){
  this_layer <- copy(stacked_layers[[idx]])
  this_layer[this_layer<0.01] <- 0
  this_layer[this_layer>0] <- 1
  return(this_layer)
})
has_transmission <- stack(has_transmission)
# levelplot(has_transmission, par.settings=rasterTheme(pal), margin=F, # at=breaks,
#           xlab=NULL, ylab=NULL, scales=list(draw=F))

par <- has_transmission * pop
names(par) <-  c(base_label, interventions)
par_millions <- cellStats(par, sum)/1000000

pdf(file.path(out_dir, paste0("atsb_par_", continent, ".pdf")), width=12, height=6)
print(levelplot(par[[1]] + 1, par.settings=rasterTheme(pal), zscaleLog=T, margin=F, # at=breaks,
                xlab=NULL, ylab=NULL, scales=list(draw=F)) +
        latticeExtra::layer(sp.polygons(africa_shp)))

print(levelplot(par[[2:5]] + 1, par.settings=rasterTheme(pal), zscaleLog=T, margin=F, # at=breaks,
                xlab=NULL, ylab=NULL, scales=list(draw=F))+
        latticeExtra::layer(sp.polygons(africa_shp)))
graphics.off()

