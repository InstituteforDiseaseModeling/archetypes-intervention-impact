###############################################################################################################
## 03_access_dev_TESTOLD.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt's original code to run the INLA model for access deviation (pixel-level
## deviation from national mean of net access).
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190430_replicate_sam joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data func_dir=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/03_access_dev_TESTOLD.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor/ --command 'Rscript ${CODE}'


rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190503_total_refactor/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190503_total_refactor/"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/run_on_gcloud/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}



## Load data and covariates
load(file.path(input_dir, "02_covariates.Rdata"))


# reload true directory names
if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190503_total_refactor/"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190503_total_refactor/"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/run_on_gcloud/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

output_fname <- file.path(output_dir, "03_access_deviation_TESTOLD.Rdata")

# load relevant functions
source(file.path(func_dir, "acc_deviation_functions.r"))


### Load (what kind of?) data ----------------------------------------------------------------------------#######################  


# same "data" as above?
data<-data[complete.cases(all.covs),]

# keep only rows with all entries filled
all.covs<-all.covs[complete.cases(all.covs),]
all.covs<-as.data.frame(all.covs)
covariate.names<-colnames(all.covs)

# merge data and covariates
data<-cbind(data,all.covs)
data<-data[sample(1:nrow(data),replace=F),]


# check for collinearity
subs<-covariate.names #<-precipitation
subs<-subs[-c(9,33,25,11)] 

X<-data[,subs]
X=cor(as.matrix(X))
diag(X)=0
ind<-which((abs(X)>0.7),arr.ind = T)

paste(subs[ind[,1]],'----',subs[ind[,2]])

STANDARD<-F
if(STANDARD){
  ms<-colMeans(data[,covariate.names])
  sds<-apply(data[,covariate.names],2,sd)
  
  for(i in 1:length(covariate.names)){
    data[,covariate.names[i]]<-(data[,covariate.names[i]]-ms[i])/(2*sds[i])
  }
}


### Run model ----------------------------------------------------------------------------#######################  

drow<-nrow(data) ## Number of data points
# load country rasters
cn<-raster(file.path(joint_dir, 'african_cn5km_2013_no_disputes.tif'))
NAvalue(cn)<--9999
POPULATIONS<-read.csv(file.path(joint_dir, 'country_table_populations.csv')) # load table to match gaul codes to country names

# adjust timing of data? this feels important****
data$yearqtr[data$yearqtr>=2015]=2014.75

# calculate access deviation for data points
data$accdev<-emplogit2(data$P,data$N)-emplogit(data$Amean,1000)
theta<-optimise(IHS.loglik, lower=0.001, upper=50, x=data$accdev, maximum=TRUE) # 2.0407
theta<-theta$maximum
data$accdev=IHS(data$accdev,theta)

# run fixed-effects linear regression of access deviation on all covariates
formula<- as.formula(paste(paste("accdev ~ "),paste(covariate.names,collapse='+'),sep=""))
l<-lm(formula,data=data)

# transform fitted values + data ... is this finding the gap? smoothing?
# what is the purpose of running this model at all?
lpacc<-plogis(emplogit(data$Amean,1000)+l$fitted.values)
#plot(data$P/data$N,lpacc)
cor(data$P/data$N,lpacc)
mean(abs((data$P/data$N-lpacc)))

# initialize inla
INLA:::inla.dynload.workaround() 

# transform data from latlong to cartesian coordinates
xyz<-ll.to.xyz(cbind(data$lon,data$lat))
data<-cbind(data,xyz)

# make a copy of the data, remove duplicates
data.est<-data
un<-paste(data.est$x,data.est$y,data.est$z,sep=':')
dup<-!duplicated(un)

# generate spatial mesh
mesh = inla.mesh.2d(loc=cbind(data.est[dup,'x'],data.est[dup,'y'],data.est[dup,'z']),
                    cutoff=0.006,
                    min.angle=c(25,25),
                    max.edge=c(0.06,500) )
print(paste('this is the number of mesh verticies',mesh$n))

# generate matern model from mesh
spde=inla.spde2.matern(mesh,alpha=2)

# generate temporal mesh
mesh1d=inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2)

# prep data for model fitting
data.cov<-data.est[,covariate.names]
est.cov<-as.list(data.cov)
est.cov$year<-data.est$yearqtr
est.cov$gaul<-data$gaul
est.cov$empmean<-data$empmean

# generate observation matrix (?)
A.est =
  inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'yearqtr'],group.mesh=mesh1d)
field.indices = inla.spde.make.index("field", n.spde=mesh$n,n.group=mesh1d$m)

# what is this?
stack.est = inla.stack(data=list(response=data.est$accdev),
                       A=list(A.est,1),
                       effects=
                         list(c(field.indices,
                                list(Intercept=1)),
                              c(est.cov)),
                       tag="est", remove.unused=TRUE)
stack.est<-inla.stack(stack.est)

formula1<- as.formula(paste(
  paste("response ~ -1 + Intercept  + "),
  paste("f(field, model=spde,group=field.group, control.group=list(model='ar1')) + ",sep=""),
  paste(covariate.names,collapse='+'),
  sep=""))


#-- Call INLA and get results --#
mod.pred =   inla(formula1,
                  data=inla.stack.data(stack.est),
                  family=c("gaussian"),
                  control.predictor=list(A=inla.stack.A(stack.est), compute=TRUE,quantiles=NULL),
                  control.compute=list(cpo=TRUE,waic=TRUE),
                  keep=FALSE, verbose=TRUE,
                  control.inla= list(strategy = 'gaussian',
                                     int.strategy='ccd', # close composite design ?
                                     verbose=TRUE,
                                     step.factor=1,
                                     stupid.search=FALSE)
)

print(paste("Saving outputs to", output_fname))
save.image(output_fname)