
library(data.table)
library(raster)
library(rasterVis)

rm(list=ls())

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results/apply_lookup")
setwd(func_dir)
source("pr_to_r0.r")
out_dir <- "output"

### Lookup table -----------------------------------------------------

# read in lookup table
lut <- fread("lookup_full_interactions_v4.csv")

# what interventions from the lookup table do you want to run?
interventions <- c("ITN 0.8; IRS 0.2;",
                   "ITN 0.8; IRS 0.8; ACT 0.8;") 

# keep only the desired interventions
lut <- lut[Intervention %in% interventions]


### Baseline rasters -----------------------------------------------------

# Make a list of any rasters for which you don't want to apply the lookup table
# This MUST include the label of the "baseline_for_lookup" raster, but can also include others
baseline_for_lookup <- "Megatrends Base 2016"
base_labels <- list("True PfPR 2017"="PfPR_rmean_Global_admin0_2017.tif",
                    "Megatrends Base 2016"="actual_ssp2_base2016_2050.tif")

baseline_pr <- raster(base_labels[[baseline_for_lookup]])

### Calculate reproductive number, or just PR?  -----------------------------------------------------

# set 'repro' to true if you want a reproductive number rather than a pfpr estimate
repro <- F
if (repro){
  print("getting repro spline")
  repro_spline <- R2spline()
}

### Desired Geography  -----------------------------------------------------

# area for which you want to predict. currently responds only to "africa" or "global"
region <- "global"

# read in cluster rasters, assign values from LUT
cluster_map <- raster(paste0(region, "_clusters_v4.tif"))

# map of cluster numbers to site names in the lut
cluster_list <-  list("aba"=1,
                      "kasama"=3,
                      "djibo"=4,
                      "kananga"=2,
                      "moine"=6,
                      "gode"=5)
if (region=="global"){
  cluster_list <- c(cluster_list, list("bajonapo"=7, "karen"=8))
}


### Functions  #####----------------------------------------------------------------------------------------------------------------------------------
get_splinefunction<-function(X,Y){
  
  plot(X,Y,type="b")
  
  splObj<-smooth.spline(X, Y)
  
  Xpred<-seq(0,1,length=100)
  Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
  # lines(Xpred,Ypred,col=2)
  return(splObj)
}

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



### Prep: Format rasters, apply R0 to baselines  #####----------------------------------------------------------------------------------------------------------------------------------

# make sure cluster map and baseline rasters are aligned
cluster_map <- crop(cluster_map, baseline_pr)
baseline_pr <- crop(baseline_pr, cluster_map)
baseline_pr <- raster::mask(baseline_pr, cluster_map)

# initialize list
raster_list <- list()
raster_idx <- 1


# loop through baseline rasters, cropping and (if necessary) converting to Rc
for (label in names(base_labels)){
  print(label)
  this_base_raster <- raster(base_labels[[label]])
  this_base_raster <- extend(crop(this_base_raster, baseline_pr), baseline_pr)
  
  if (repro==T){
    print("finding rc")
    this_base_rc <-apply_spline_to_raster(repro_spline, this_base_raster, return_r0=T)
    raster_list[[raster_idx]] <- this_base_rc
  }else{
    raster_list[[raster_idx]] <- this_base_raster
  }
  raster_idx <- raster_idx + 1
}

### Main loop: apply lookup table  #####----------------------------------------------------------------------------------------------------------------------------------

for (intervention in interventions){
  print(intervention)
  
  # loop through each cluster, apply the appropriate spline, convert to reproductive number if desired
  final_prs <- lapply(names(cluster_list), function(site_name){
    this_mask <- cluster_map==cluster_list[[site_name]]
    pr_masked <- raster::mask(baseline_pr, this_mask, maskvalue=0)
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
  # this_pr_final <- mask(this_pr_final, cluster_map, maskvalue=FALSE)
  raster_list[[raster_idx]] <- this_pr_final
  raster_idx <- raster_idx + 1 
}

stacked_layers <- stack(raster_list)
names(stacked_layers) <- c(names(base_labels), interventions)

### Plotting and saving  #####----------------------------------------------------------------------------------------------------------------------------------

if (repro==T){
  breaks <- c(0,1,seq(1.5, 5, length.out=25), seq(5.1, 80, length.out=4))
  pal = c("#e0e0e0", terrain.colors(31)[1:30])
  pdf(file.path(out_dir, paste0("RC_", region, ".pdf")), width=12, height=6)
  print(levelplot(stacked_layers, par.settings=rasterTheme(pal), at=breaks, xlab=NULL, ylab=NULL, scales=list(draw=F)))
  graphics.off()
  writeRaster(stacked_layers, options="INTERLEAVE=BAND", bylayer=T, suffix="names", filename=file.path(out_dir, paste0("RC_",region, ".tif")), overwrite=T)
}else{
  pal <- c("#e0e0e0", rev(brewer.pal(11, "Spectral")))
  breaks <- c(0, seq(0.01, 1, length.out=11))
  pdf(file.path(out_dir, paste0("pfpr_", region, ".pdf")), width=12, height=6)
  print(levelplot(stacked_layers, par.settings=rasterTheme(pal), at=breaks, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F)))
  graphics.off()
  writeRaster(stacked_layers, options="INTERLEAVE=BAND", bylayer=T, suffix="names", filename=file.path(out_dir, paste0("pfpr_",region, ".tif")), overwrite=T)
}


