###############################################################################################################
## 03_model_access_dev.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt's original code to run the INLA model for access deviation (pixel-level
## deviation from national mean of net access).
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/02_prep_covariates.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor/ --command 'Rscript ${CODE}'


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
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/joint_data"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

# load relevant functions
source(file.path(func_dir, "03_access_dev_functions.r"))
# go through script, remake functions as you go
# source(file.path(func_dir, "acc_deviation_functions.r"))
# source(file.path(func_dir, "algorithm.V1.R"))
# source(file.path(func_dir, "INLAFunctions.R"))

output_fname <- file.path(output_dir, "03_access_deviation.Rdata")


## Load data 
data <- fread(file.path(input_dir, "02_covariates.csv"))

## Prep for modeling ## ---------------------------------------------------------

# drop NAs
dropped_rows <- nrow(data) - nrow(data[complete.cases(data)])
print(paste("dropping", dropped_rows, "rows of data due to null values in covariates!"))
data <- data[complete.cases(data)]

# shuffle row order (why?)
data <- data[sample(1:nrow(data),replace=F),]

# check for collinearity
cov_names <- names(data)[(which(names(data)=="row_id")+1):length(names(data))]
cov_data <- data[, cov_names, with=F]
collin <- cor(as.matrix(cov_data))
diag(collin) <- 0
high_collin <- which(abs(collin)>0.7, arr.ind=T)

if (nrow(high_collin)>0){
  warning("Collinear covariates identified!")
  print(high_collin)
}


## Run model ##-------------------------------------------------------------

# adjust timing of data? this feels important****
data[, yearqtr:=min(yearqtr, 2014.75)]


# calculate access deviation for data points, 
# transform via empirical logit and inverse hyperbolic sine
data[, emp_access_dev:= emplogit2(P, N) - emplogit(Amean, 1000)]

theta<-optimise(IHS.loglik, lower=0.001, upper=50, x=data$emp_access_dev, maximum=TRUE) # 2.0407
theta<-theta$maximum

data[, ihs_emp_access_dev:=IHS(emp_access_dev, theta)] 


# transform data from latlong to cartesian coordinates
xyz<-ll.to.xyz(data[, list(row_id, longitude=lon, latitude=lat)])
data <- merge(data, xyz, by="row_id", all=T)

# initialize inla
INLA:::inla.dynload.workaround() 

# make a copy of the data for estimation
data_est <- copy(data)

# generate spatial mesh using unique xyz values 
# TODO: is this all xyz is used for?
mesh = inla.mesh.2d(loc= unique(data_est[, list(x,y,z)]),
                    cutoff=0.006,
                    min.angle=c(25,25),
                    max.edge=c(0.06,500) )
print(paste("Mesh constructed:", mesh$n, "vertices"))

# generate spde matern model from mesh
spde_matern =inla.spde2.matern(mesh,alpha=2) # formerly spde

# generate temporal mesh
temporal_mesh=inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2) # formerly mesh1d


### un-refactored from here down

# prep data for model fitting
data.cov<-data.est[,covariate.names]
est.cov<-as.list(data.cov)
est.cov$year<-data.est$yearqtr
est.cov$gaul<-data$gaul
est.cov$empmean<-data$empmean

# generate observation matrix (?)
A.est =
  inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'yearqtr'],group.mesh=temporal_mesh)
field.indices = inla.spde.make.index("field", n.spde=mesh$n,n.group=temporal_mesh$m)

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
inla.ks.plot(mod.pred$cpo$pit, punif)




