#################################################################################################
get_splinefunction_from_LUT<-function(X,Y){
    
    plot(X,Y,type="b")
    
    #lmObj<- lm(Y~poly(X,ORDER,raw=TRUE))
    splObj<-smooth.spline(X, Y)

    Xpred<-seq(0,1,length=100)
    Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
    lines(Xpred,Ypred,col=2)
    return(splObj)
}
##################################################################################################
apply_LUT_to_raster<-function(splObj,inRaster){

    inMat<-as.matrix(inRaster)
    inVec<-as.vector(inMat)
    NAid<-is.na(inVec)
    inVecNoNA<-inVec[!NAid]

    outVecNoNA<-predict(splObj,inVecNoNA)$y
    outVec<-rep(NA,length(inVec))
    outVec[!NAid]<-outVecNoNA
    outMat<-matrix(outVec,nrow=nrow(inMat))
    outRaster<-inRaster
    outRaster[]<-outMat

    return(outRaster)
}
##################################################################################################
library(raster)

##
## MEG_may2018_pilot_NO INTERVENTIONS BASELINE_with_LUTv2 
##

# import LUT
LUT<-read.csv("C:/PROJECTS/70_Megatrends/analysis/pilot_SAGme_june2018/LUTs/lookup_table_v2.csv")
names(LUT)

for(YEAR in c(2016,2020,2030,2040,2050)){
    inRasterPATH<-paste0("C:/PROJECTS/70_Megatrends/analysis/pilot_SAGme_june2018/output/PfPR/MEG_may2018_pilot_NO INTERVENTIONS BASELINE/MEGAT_",YEAR,"_6.tif")
    outRasterPATH<-paste0("C:/PROJECTS/70_Megatrends/analysis/pilot_SAGme_june2018/output/PfPR/MEG_may2018_pilot_NO INTERVENTIONS BASELINE_with_LUTv2/MEGAT_",YEAR,"_6_LUTv2.tif")
    inRaster<-raster(inRasterPATH)
    splObj<-get_splinefunction_from_LUT(LUT$initial,LUT$final)
    outRaster<-apply_LUT_to_raster(splObj,inRaster)
    writeRaster(outRaster,outRasterPATH)
}



