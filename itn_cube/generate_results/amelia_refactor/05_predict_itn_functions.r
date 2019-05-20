###############################################################################################################
## 05_predict_itn_functions.r
## Amelia Bertozzi-Villa
## May 2019
## 
## Functions to accompany 05_predict_itn.r-- prepping survey data to go into the ITN cube model. 
##############################################################################################################


ll.to.xyz.old<-function(ll){
  
  if(is.null(colnames(ll))){
    colnames(ll)<-c('longitude','latitude')
  }
  if(colnames(ll)[1]=='x' & colnames(ll)[2]=='y'){
    colnames(ll)<-c('longitude','latitude')
  }
  if(colnames(ll)[1]=='lon' & colnames(ll)[2]=='lat'){
    colnames(ll)<-c('longitude','latitude')
  }
  ll[,'longitude']<-ll[,'longitude']*(pi/180)
  ll[,'latitude']<-ll[,'latitude']*(pi/180)
  x = cos(ll[,'latitude']) * cos(ll[,'longitude'])
  y = cos(ll[,'latitude']) * sin(ll[,'longitude'])
  z = sin(ll[,'latitude'])
  return(cbind(x,y,z))
}