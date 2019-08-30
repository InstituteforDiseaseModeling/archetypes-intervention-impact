
library(dichromat)
library(MapSuite)
library(raster)





### Spline application Functions   #####----------------------------------------------------------------------------------------------------------------------------------
get_splinefunction<-function(X,Y){
  
  # plot(X,Y,type="b")
  
  splObj<-smooth.spline(X, Y)
  
  Xpred<-seq(0,1,length=100)
  Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
  lines(Xpred,Ypred,col=2)
  return(splObj)
}

apply_spline_to_raster<-function(splObj,inRaster, return_r0=F){ # todo: there must be a better way
  
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


apply_lookup <- function(intervention, baseline_pr, lut, cluster_map, bounding_pr){
  final_prs <- lapply(unique(lut$Site_Name), function(site_name){
    this_mask <- cluster_map==unique(lut[Site_Name==site_name]$cluster_id)
    pr_masked <- raster::mask(baseline_pr, this_mask, maskvalue=0)
    this_lut <- lut[Site_Name==site_name & Intervention==intervention, list(mean_initial, mean_final)]
    this_spline <- get_splinefunction(this_lut$mean_initial,this_lut$mean_final)
    this_pr <- apply_spline_to_raster(this_spline, pr_masked, return_r0 = F)
    
    return(this_pr)
    
  })
  
  # merge back into a single raster
  this_pr_final <- do.call(merge, final_prs)
  
  # bound if needed
  this_pr_final <- min(stack(this_pr_final, bounding_pr)) 
  names(this_pr_final) <- intervention
  
  return(this_pr_final)
}

### Plotting and color functions   #####----------------------------------------------------------------------------------------------------------------------------------

generate_ramp <- function(colors, vals, n.out=10){
  rampfunc <- colorRampPalette(colors)
  pal <- rampfunc(n.out)
  breaks <- seq(vals[1], vals[2], length.out=n.out)
  return(list(pal, breaks))
}

generate_full_pal <- function(){
  
  map_pal <-c("#E6E6E6", "#A3A3A3", "#5A17FD", "#51C8FE", "#AEFEAB", "#FBD817",  "#FE3919") 
  map_breaks <- c(0, 0.005, 0.01, 0.25, 0.5, 0.75, 1)
  
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
  
  return(list(pal=full_pal, breaks=full_breaks))
}

