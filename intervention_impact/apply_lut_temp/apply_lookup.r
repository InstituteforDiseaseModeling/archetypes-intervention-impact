
library(data.table)
library(raster)
library(rasterVis)

rm(list=ls())

source("pr_to_r0.r")

main_dir <- "/Users/bertozzivill/Desktop/apply_lut_temp"
interventions <- c("ITN 0.4; ", "ITN 0.4; IRS 0.4; ", "ITN 0.4; ACT 0.6; ") # todo: remove final space from intervention name
repro <- T

#################################################################################################
get_splinefunction_from_LUT<-function(X,Y){
  
  # plot(X,Y,type="b")
  
  #lmObj<- lm(Y~poly(X,ORDER,raw=TRUE))
  splObj<-smooth.spline(X, Y)
  # 
  # Xpred<-seq(0,1,length=100)
  # Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
  # lines(Xpred,Ypred,col=2)
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
  outVecNoNA <- pmin(outVecNoNA, inVecNoNA)
  outVecNoNA[outVecNoNA<0] <- 0
  
  if (return_r0){
    outVecNoNA <- unlist(lapply(outVecNoNA, PR2R))
  }
  
  outVec<-rep(NA,length(inVec))
  outVec[!NAid]<-outVecNoNA
  outMat<-matrix(outVec,nrow=nrow(inMat))
  outRaster<-inRaster
  outRaster[]<-outMat
  
  return(outRaster)
}
##################################################################################################

lut <- fread(file.path(main_dir, "lookup_full_interactions.csv"))
lut <- lut[Intervention %in% interventions]

# view plots-- Karen's is weird at high prev but we shouldn't need it up there.
# for (site_name in unique(lut$Site_Name)){
#   test <- unique(lut[Site_Name==site_name, list(mean_initial, mean_final)])
#   get_splinefunction_from_LUT(test$mean_initial,test$mean_final)
# }

# read in rasters: PR and clustering
pr_orig <- raster(file.path(main_dir, "global_PfPR.2016.tif"))
masks <- raster(file.path(main_dir, "MAP_Regions_Pf_5k.tif"))
cluster_map <- raster(file.path(main_dir, "clusters_joint_6.tif"))

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

pdf("rc_applied.pdf", width=9, height=6)

for (continent in names(mask_values)){
  
  print(continent)
  
  if (continent=="africa"){
      continent_mask <- copy(cluster_map)
      maskvalue=0
      
      final_mask <- masks == mask_values[[continent]]
      final_mask <- trim(final_mask, values=F)
      
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
    null_spline <- get_splinefunction_from_LUT(seq(0,1,length=100), rep(1, 100))
    r0 <-apply_LUT_to_raster(null_spline, this_pr_orig, return_r0 = T)
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
      this_pr <- apply_LUT_to_raster(this_spline, pr_masked, return_r0 = T)
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
    names(stacked_layers) <- c("R0 2016", interventions)
    breaks <- c(seq(0,1,length.out=50), seq(1,10, length.out = 26)[2:26], seq(10, 45, length.out = 26)[2:26])
    print(levelplot(stacked_layers, par.settings=rasterTheme(blue2red(30)), at=breaks))
  }else{
    names(stacked_layers) <- c("PfPR 2016", interventions)
    print(levelplot(stacked_layers, par.settings=rasterTheme(brewer.pal(7, "BuPu"))))
  }
  
}

graphics.off()

# map_plot <- levelplot(cluster_raster, att="ID", col.regions=brewer.pal(nclust, palette),
#                       xlab=NULL, ylab=NULL, scales=list(draw=F))
#                       

