###############################################################################################################
## scratch code for raking to stock and flow-- add to 05_predict_itn when done 
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190521_replicate_prediction func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data/ --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/05_predict_itn.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190521_replicate_prediction/05_predictions --command 'Rscript ${CODE}'

rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo", "VGAM", "raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190521_replicate_prediction//"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190521_replicate_prediction//"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
  cov_dir <- "/Volumes/GoogleDrive/Team Drives/cubes/5km incomplete/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
  cov_dir <- Sys.getenv("cov_dir")
}

# directories
joint_dir <- file.path(joint_dir, "For_Prediction")
set.seed(212)
prediction_years <- 2000:2016

## -------------

print("loading indicators") # todo: remind myself what these actually are?
quarterly_inds <- fread(file.path(joint_dir, 'indicators_access_qtr_new.csv'))
quarterly_inds <- quarterly_inds[2:nrow(quarterly_inds)]
names(quarterly_inds) <- c("country", seq(2000, 2018, 0.25))
quarterly_inds <- melt(quarterly_inds, id.vars = "country", variable.name="year")
quarterly_inds[, year:=floor(as.numeric(as.character(year)))]

annual_inds <- quarterly_inds[year<=max(prediction_years), list(value=mean(value)), by=list(country, year)]


### OLD: 

print("loading indicators")
indicators<-read.csv(file.path(joint_dir, 'indicators_access_qtr_new.csv'))


ind<-matrix(nrow=nrow(indicators),ncol=17)
for(i in 1:nrow(indicators)){
  ind[i,]<-tapply(as.numeric(indicators[i,2:ncol(indicators)]),c(rep(1:(72/4), each = 4),18),mean)[1:17]
}
rownames(ind)<-indicators[,1]
colnames(ind)<-seq(2000,2016,by=1)
indicators<-ind

# todo: does this POPULATIONS pull a different file than the one earlier in the script?
print("loading populations")
POPULATIONS<-read.csv(file.path(joint_dir, '../country_table_populations.csv')) # load table to match gaul codes to country names
names<-as.character(rownames(ind))
#for(i in 1:nrow(indicators)){
#	names[i]=as.character(POPULATIONS[as.character(POPULATIONS$NAME)==names[i],'COUNTRY_ID']) # get 3 letter country codes
#}
print("generating indicator matrix")
rownames(ind)<-as.character(names)
threshold=0.02
indmat<-matrix(data=0,nrow=nrow(indicators),ncol=ncol(indicators))
for(i in 1:nrow(indicators)){
  for(j in 1:ncol(indicators)){
    if(indicators[i,j]<threshold){
      indmat[i,j]=1
    }
    if(indicators[i,j]>threshold){
      break
    }
  }
}
colnames(indmat)<-seq(2000,2016,by=1)
rownames(indmat)<-names

times<-seq(2000,2016,by=1)

fnames_for_stack <- paste(out_dir, '/ITN_',times,'.USE.tif',sep="")
print("contents of stack directory:")
print(out_dir)
list.files(out_dir)
print("loading raster stack")
print(fnames_for_stack)

st<-stack(paste(out_dir, '/ITN_',times,'.USE.tif',sep=""))
print("plotting original use stack")

pal <- brewer.pal(8, "RdYlGn")
breaks <- seq(0, 1, length.out = length(pal)+1)

pdf(file.path(out_dir, "samdir_zdir_compare.pdf"), width=7, height=10)
p <- levelplot(st,
               par.settings=rasterTheme(region=pal), at=breaks,
               xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F,
               main="Initial Use Rasters"
)
print(p)


stv<-getValues(st)
cn<-raster(file.path(joint_dir, '../african_cn5km_2013_no_disputes.tif'))
cnv<-getValues(cn)
for(i in 1:nrow(indmat)){
  gaul<-POPULATIONS[as.character(POPULATIONS$NAME)==rownames(ind)[i],'GAUL_CODE']
  tmp<-cn
  wh<-colnames(indmat)[indmat[i,]==1]
  wh<-as.numeric(wh)-1999
  stv[cnv==gaul,wh]<-0
}

for(i in 1:nlayers(st)){
  st[[i]]<-setValues(st[[i]],stv[,i])
}

p <- levelplot(st,
               par.settings=rasterTheme(region=pal), at=breaks,
               xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F,
               main="Use Rasters for Z Drive"
)
print(p)
graphics.off()

print("saving for z dir rasters")

for(i in 1:nlayers(st)){
  writeRaster(st[[i]],file.path(out_dir, paste("for_z_", times[i],'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
  
}


