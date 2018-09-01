
library(data.table)
library(raster)
library(rasterVis)
library(colorRamps)
library(MapSuite)

rm(list=ls())

source("pr_to_r0.r")


interventions <- c("ACT 0.2; ", "IRS 0.4; ACT 0.2; " , "ITN 0.4; IRS 0.4; ACT 0.2; ") # todo: remove final space from intervention name
repro <- T
if (repro){
  print("getting repro spline")
  repro_spline <- R2spline()
}

#################################################################################################
get_splinefunction_from_LUT<-function(X,Y){
  
  plot(X,Y,type="b")
  
  splObj<-smooth.spline(X, Y)
 
  Xpred<-seq(0,1,length=100)
  Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
  lines(Xpred,Ypred,col=2)
  return(splObj)
}
##################################################################################################

##################################################################################################
apply_LUT_to_raster<-function(splObj,inRaster, return_r0=F){
  
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

lut <- fread(file.path("..", "shiny", "data", "lookup_full_interactions_v2.csv"))
lut <- lut[Intervention %in% interventions]

# view plots-- Karen's is weird at high prev but we shouldn't need it up there.
# for (site_name in unique(lut$Site_Name)){
#   test <- unique(lut[Site_Name==site_name, list(mean_initial, mean_final)])
#   get_splinefunction_from_LUT(test$mean_initial,test$mean_final)
# }

# read in rasters: PR and clustering
pr_orig <- raster(file.path("MODEL43.2015.PR.ALL.rmean.tif"))

masks <- raster("MAP_Regions_Pf_5k.tif")
cluster_map <- raster("clusters_joint_6.tif")
# africa counterfactual is smaller than cluster map
cluster_map <- crop(cluster_map, pr_orig)
pr_orig <- crop(pr_orig, cluster_map)

mask_values <- list("africa"=1,
                    "americas"=2,
                    "asia"=3)
mask_values <- list("africa"=1)

cluster_list <- list("africa"= list("aba"=1,
                                "martae"=2,
                                "djibo"=3,
                                "kananga"=4,
                                "moine"=5,
                                "gode"=6),
                 "americas"=list("bajonapo"=1),
                 "asia" = list("karen"=1))



for (continent in names(mask_values)){
  
  print(continent)
  
  if (continent=="africa"){
      continent_mask <- copy(cluster_map)
      maskvalue=0
      
      final_mask <- masks == mask_values[[continent]]
      final_mask <- trim(final_mask, values=F)
      
      # again, crop mask to smaller africa raster
      final_mask <- crop(final_mask, pr_orig)
      
    }else{
      continent_mask <- masks == mask_values[[continent]]
      continent_mask <- trim(continent_mask, values=F)
      maskvalue=F
      
      final_mask <- copy(continent_mask)
    } 
  
  this_pr_orig <- crop(pr_orig, final_mask)
  this_pr_orig <- mask(this_pr_orig, final_mask, maskvalue=F)
  
  if (repro==T){
    print("finding r0 for raster")
    r0 <-apply_LUT_to_raster(repro_spline, this_pr_orig, return_r0=T)
    raster_list <- list(r0)
  }else{
    raster_list <- list(this_pr_orig)
  }
  
  idx <- 2
  
  for (intervention in interventions){
    print(intervention)
    
    final_prs <- lapply(names(cluster_list[[continent]]), function(site_name){
      this_mask <- continent_mask==cluster_list[[continent]][[site_name]]
      pr_masked <- mask(this_pr_orig, this_mask, maskvalue=maskvalue)
      this_lut <- unique(lut[Site_Name==site_name & Intervention==intervention, list(mean_initial, mean_final)])
      this_spline <- get_splinefunction_from_LUT(this_lut$mean_initial,this_lut$mean_final)
      this_pr <- apply_LUT_to_raster(this_spline, pr_masked, return_r0 = F)
      this_reptro <- apply_LUT_to_raster(repro_spline, this_pr, return_r0 = T)
    })
    
    if (length(final_prs)>1){
      this_pr_final <- do.call(merge, final_prs)
      this_pr_final <- mask(this_pr_final, final_mask, maskvalue=FALSE)
      
    }else{
      this_pr_final <- final_prs[[1]]
    }
    
    raster_list[[idx]] <- this_pr_final
    idx <- idx + 1 
  }
  
  
  stacked_layers <- stack(raster_list)
  if (repro==T){
    names(stacked_layers) <- c("R0 2015", interventions)
    breaks <- c(seq(-0.001, 5, length.out=50), seq(5, 10, length.out = 46)[2:46], seq(10, 75, length.out=6)[2:6])
    pdf("rc_applied.pdf", width=9, height=6)
    print(levelplot(stacked_layers, par.settings=rasterTheme(region=wpal("sky")), zscaleLog=T, xlab=NULL, ylab=NULL, scales=list(draw=F)))
    graphics.off()
  }else{
    names(stacked_layers) <- c("PfPR 2015", interventions)
    print(levelplot(stacked_layers, par.settings=rasterTheme(brewer.pal(7, "BuPu"))))
  }
  
}

writeRaster(stacked_layers, options="INTERLEAVE=BAND", "repro_numbers.tif", overwrite=T)


# map_plot <- levelplot(cluster_raster, att="ID", col.regions=brewer.pal(nclust, palette),
#                       xlab=NULL, ylab=NULL, scales=list(draw=F))
#                       

