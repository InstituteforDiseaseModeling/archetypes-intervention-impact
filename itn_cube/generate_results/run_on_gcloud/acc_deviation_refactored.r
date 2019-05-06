

# unused covariates:
# '/home/drive/cubes/5km/LST_delta/mean/' #
# '/home/drive/cubes/5km/TCB/mean/' #
# '/home/drive/cubes/5km/Topographic/Africa_FA_90m.mean.tif' #
# '/home/drive/cubes/5km/worldclim/prec57a0.tif'


# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/cubes_5km joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data func_dir=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud --input database_fname=gs://map_data_z/users/amelia/itn_cube/results/20190430_replicate_sam/01_database.csv CODE=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud/acc_deviation_refactored.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190430_replicate_sam/ --command 'Rscript ${CODE}'



rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/joint_data"
  database_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/output/ITN_final_clean_access_9Feb2019.csv"
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/input"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/access_deviation"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  database_fname <- Sys.getenv("database_fname") # location of output file from generate_database_refactored.r
  input_dir <- Sys.getenv("input_dir") # here, location of covariate data 
  output_dir <- Sys.getenv("output_dir") 
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

# load relevant functions
source(file.path(func_dir, "acc_deviation_functions.r"))
source(file.path(func_dir, "algorithm.V1.R"))
source(file.path(func_dir, "INLAFunctions.R"))

set.seed(212)

covariate_fname <- file.path(output_dir, "02_covariates_TESTREPLICATION_TAKETWO.Rdata")
output_fname <- file.path(output_dir, '03_access_deviation_TESTREPLICATION_TAKETWO.Rdata')

# Load data from create_database.r  ------------------------------------------------------------
data<-read.csv(database_fname)
data<-data[!is.na(data$year),]

### Load covariates  ----------------------------------------------------------------------------#######################  

# Load dynamic covariates  ------------------------------------------------------------

dynamic_covnames <- c("LST_day/mean", "LST_night/mean", "EVI/mean", "TCW/mean", "TSI/mean")



#' foldersd<-c('/home/drive/cubes/5km/LST_day/mean/',
#'             #'/home/drive/cubes/5km/LST_delta/mean/',
#'             '/home/drive/cubes/5km/LST_night/mean/',
#'             '/home/drive/cubes/5km/EVI/mean/',
#'             #'/home/drive/cubes/5km/TCB/mean/',
#'             '/home/drive/cubes/5km/TCW/mean/',
#'             '/home/drive/cubes/5km/TSI/mean/')

# find raster cell ids from access/use data 
dates<-as.Date(as.yearmon(data$year))
split.dates<-str_split_fixed(dates, "-", 3)
uniques<-unique(split.dates)
l=length(dynamic_covnames)
cell<-data$cellnumber

# extract covariates for unique year-months	

ncores <- detectCores()
print(paste("--> Machine has", ncores, "cores available"))
registerDoParallel(ncores-2)
print("Extracting year-month covariates")
covs.list.dyn<-foreach(i=1:nrow(uniques)) %dopar% { # loop through unique names
  print(paste("unique date", i, "of", nrow(uniques)))
  wh<-split.dates[,1]==uniques[i,1] & split.dates[,2]==uniques[i,2]
  print(uniques[i,])
  un.cells<-cell[wh]
  tmp<-matrix(NA,nrow=length(un.cells),ncol=l)
  for(j in 1:l){
    year<-as.numeric(uniques[i,1])
    cov_dir <- file.path(input_dir, dynamic_covnames[j])
    
    while(!file.exists(paste0(cov_dir, "/", year,'.',uniques[i,2],'.mean.tif'))){	# some years have no data yet
      year=year-1
    }
    
    this_file_path <- paste0(cov_dir, "/", year,'.',uniques[i,2],'.mean.tif')
    
    r=raster(this_file_path)
    NAvalue(r)=-9999
    tmp[,j]<-r[un.cells]
  }
  return(tmp)
}

# assign these to the correct cells in data
dynamic.covs<-matrix(NA,nrow=nrow(data),ncol=length(dynamic_covnames))
for(i in 1:nrow(uniques)){
  wh<-split.dates[,1]==uniques[i,1] & split.dates[,2]==uniques[i,2]
  dynamic.covs[wh,]<-covs.list.dyn[[i]]
}
#colnames(dynamic.covs)<-c('lst_day','lst_delta','lst_night','evy','tcb','tcw','tsi')
colnames(dynamic.covs)<-c('lst_day','lst_night','evy','tcw','tsi')
print("year-month covariates extracted successfully")

# Load year-only (??) covariates  ------------------------------------------------------------
# foldery<-c('/home/drive/cubes/5km/IGBP_Landcover/Fraction/','/home/drive/cubes/5km/AfriPop/')
year_covnames <- c("IGBP_Landcover/Fraction", "AfriPop")

uniques<-unique(split.dates[,1])

l=17+1 # 17 fraction classes and 1 afripop
cell<-data$cellnumber

#compute covariates for unique years
print("Extracting annual covariates")
registerDoParallel(ncores-2)
covs.list.year<-foreach(i=1:length(uniques)) %dopar% { # loop through unique names
  wh<-split.dates[,1]==uniques[i] 
  un.cells<-cell[wh]
  tmp<-matrix(NA,nrow=length(un.cells),ncol=l)
  for(j in 1:(l-1)){
    # no land cover for 2013
    if(uniques[i]>2012){
      r=raster(paste0(file.path(input_dir, "IGBP_Landcover/Fraction"), "/",'2012','.fraction.class.',j-1,'.tif'))
    }else{
      r=raster(paste0(file.path(input_dir, "IGBP_Landcover/Fraction"), "/",uniques[i],'.fraction.class.',j-1,'.tif'))
    }		
    NAvalue(r)=-9999
    tmp[,j]<-r[un.cells]
  }
  if(uniques[i]>2015){
    r=raster(paste0(file.path(input_dir, "AfriPop"), "/",'2015','.total.population.tif'))
  }else{
    r=raster(paste0(file.path(input_dir, "AfriPop"), "/",uniques[i],'.total.population.tif'))
  }	
  NAvalue(r)=-9999
  tmp[,l]<-r[un.cells]
  return(tmp)
}

#assign these to the correct cells in data
yearonly.covs<-matrix(NA,nrow=nrow(data),ncol=l)
for(i in 1:length(uniques)){
  wh<-split.dates[,1]==uniques[i] 
  yearonly.covs[wh,]<-covs.list.year[[i]]
}
colnames(yearonly.covs)<-c(paste0('landcover',0:16),'populatopn')

yearonly.covs<-yearonly.covs[,-14] # remove landcover 13  - Urban and built-up   - for collinearity with population
print("annual covs extracted successfully")

# Load static covariates  ------------------------------------------------------------


print("Extracting Static covariates")
static_fnames <- c("Topographic/Africa_TMI_90m.mean.tif",
                   "Topographic/Africa_SRTM_90m.mean.tif",
                   "Topographic/Africa_slope_90m.mean.tif",
                   "Seasonality/pf_seasonality.tif",
                   "Poverty/PET_1950-2000_5km.mean.tif",
                   "Poverty/AI_1950-2000_5km.mean.tif",
                   "Poverty/accessibility_50k_5km.mean.tif",
                   "Poverty/viirs_nighttime_5km.mean.tif")

static_covnames <- file.path(input_dir, static_fnames)

#' folderss<-c('/home/drive/cubes/5km/Topographic/Africa_TMI_90m.mean.tif',
#'             '/home/drive/cubes/5km/Topographic/Africa_SRTM_90m.mean.tif',
#'             '/home/drive/cubes/5km/Topographic/Africa_slope_90m.mean.tif',
#'             #'/home/drive/cubes/5km/Topographic/Africa_FA_90m.mean.tif',
#'             '/home/drive/cubes/5km/Seasonality/pf_seasonality.tif',
#'             '/home/drive/cubes/5km/Poverty/PET_1950-2000_5km.mean.tif',
#'             '/home/drive/cubes/5km/Poverty/AI_1950-2000_5km.mean.tif',
#'             '/home/drive/cubes/5km/Poverty/accessibility_50k_5km.mean.tif',
#'             #'/home/drive/cubes/5km/worldclim/prec57a0.tif',
#'             '/home/drive/cubes/5km/Poverty/viirs_nighttime_5km.mean.tif')
#' #'/home/drive/cubes/5km/Poverty/DMSP_F18_5km.mean.tif')


st<-stack(static_covnames)
NAvalue(st)=-9999

static.covs<-st[data$cellnumber]

all.covs<-cbind(static.covs,yearonly.covs,dynamic.covs)

print("Static covariates extracted successfully")

save.image(covariate_fname)


##todo: start new script here 

### Load (what kind of?) data ----------------------------------------------------------------------------#######################  

# load saved covariates
load(covariate_fname)
## IHS: BURBRIDGE-- inverse hyperbolic sine transform


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

print("fixed effects:")
print(mod.pred$summary.fixed)


print(paste("Saving outputs to", output_fname))
save.image(output_fname)
inla.ks.plot(mod.pred$cpo$pit, punif)


# test: comment out everything south of here to see if that writes the results

### Assess model ----------------------------------------------------------------------------#######################

# compute summary statistics and information criteria (posterior fit)

# watanabe-aikake information criterion-- this code seems to be expecting a list of models from which we can 
# query waics, but the above code does not generate such a list.
# waic<-c()
# for(i in 1:4){
#   waic[i]<-models[[i]]$waic$waic
# }
# a=(waic[1]-waic[2])/waic[1]
# b=(waic[2]-waic[3])/waic[1]
# c=(waic[3]-waic[4])/waic[1]
# 
# a/(a+b+c)
# b/(a+b+c)
# c/(a+b+c)

# more model stats? what is lpgap?
index= inla.stack.index(stack.est,"est")$data
lp=mod.pred$summary.linear.predictor$mean[index]
lp=Inv.IHS(lp,theta)

# calculate fit
lpgap<-plogis(emplogit(data$Amean,1000)+lp)
plot(data$P/data$N,lpgap)
cor(data$P/data$N,lpgap)
mean(abs((data$P/data$N-lpgap)))

a=mod.pred$summary.fixed$mean
sig<-sign(mod.pred$summary.fixed$'0.025quant')==sign(mod.pred$summary.fixed$'0.975quant')
names(a)<-rownames(mod.pred$summary.fixed)
sort(a[sig])

### Cross validation ----- (looc w/o re-running?)
y=data.est$accdev

Q.space = inla.spde.precision(spde, theta=c(mod.pred$summary.hyperpar$mean[2],mod.pred$summary.hyperpar$mean[3]))
a=mod.pred$summary.hyperpar$mean[4]
Q.time =
  sparseMatrix(i=c(1:mesh1d$m, 1:(mesh1d$m-1)),
               j=c(1:mesh1d$m, 2:mesh1d$m),
               x=(c(c(1, rep((1+a^2), mesh1d$m-2), 1),
                    rep(-a, mesh1d$m-1)) /(1-a^2)),
               dims=mesh1d$m*c(1,1),
               symmetric=TRUE)

# Space-time precision: the magic
Q = kronecker(Q.time, Q.space)

#A =inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']))
A= inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'yearqtr'],group.mesh=mesh1d)

# todo: work through this
index= inla.stack.index(stack.est,"est")$data
my=mod.pred$summary.linear.predictor$mean[index]-A%*%mod.pred$summary.random$field$mean
Qe<-Diagonal(n=length(my),x=mod.pred$summary.hyperpar$mean[1])
w=t(A)%*%Qe%*%(y-my)
Qxy <- Q+t(A)%*%Qe%*%A
z<- inla.qsolve(Qxy,w,method='solve')
mu_xy <- drop(my + A%*%drop(z))
S=inla.qinv(Qxy)
V<-rowSums(A * (A %*% S))
q=diag(Qe)


### Anomaly -----

# plot... more summaries? what is mval?
mval<-rep(NA,nrow(data))
for(i in 1:nrow(data)){
  mval[i]<-y[i]-(y[i]-mu_xy[i])/(1-(q[i]*V[i]))
}
mval<-Inv.IHS(mval,theta)
mval<-plogis(emplogit(data$Amean,1000)+mval)
plot(data$P/data$N,mval)
cor(data$P/data$N,mval)
mean(abs((data$P/data$N-mval)))

plot(y,mval)
plot(data$P/data$N,plogis(emplogit(data$Amean,1000)+mval))
cor(data$P/data$N,plogis(emplogit(data$Amean,1000)+mval))
mean(abs(data$P/data$N-plogis(emplogit(data$Amean,1000)+mval)))


un=unique(data$Survey)
actual<-predicted<-c()

for(i in 1:length(un)){
  # sam: hammer this home (survey-level predictions are good)
  actual[i]<-mean(data$P[data$Survey==un[i]]/data$N[data$Survey==un[i]])
  predicted[i]<-mean(plogis(emplogit(data$Amean[data$Survey==un[i]],1000)+mval[data$Survey==un[i]]))

}

plot(predicted,actual,pch=16,col=1:length(un))
abline(0,1,col='red')
cor(predicted,actual)
mean(abs(predicted-actual))
mean((predicted-actual)^2)


spde.res2 <- inla.spde2.result(mod.pred, "field", spde)
range=spde.res2$marginals.range.nominal$range.nominal.1
plot(gc.dist(range[,1]),range[,2])
variance=spde.res2$marginals.variance.nominal$variance.nominal.1
print(paste('The mean range is',gc.dist(exp(spde.res2$summary.log.range.nominal$mean))))



