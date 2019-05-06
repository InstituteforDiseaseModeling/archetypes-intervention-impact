###############################################################################################################
## 03_access_dev.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt's original code to run the INLA model for access deviation (pixel-level
## deviation from national mean of net access).
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/03_access_dev.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190503_total_refactor/ --command 'Rscript ${CODE}'


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

set.seed(212)

output_fname <- file.path(output_dir, "03_access_deviation.Rdata")


## Load data 
data <- fread(file.path(input_dir, "02_covariates.csv"))
# row_order <- fread(file.path(input_dir, "03_row_order_FROM_20190430_replicate_sam.csv"))

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

# data_est <- merge(data_est, row_order, by="row_id", all=T)
# data_est <- data_est[order(shuffled_row_location)]

# generate spatial mesh using unique xyz values 
# TODO: is this all xyz is used for?
# TODO: ask Sam about discrepancy in mesh count. 
spatial_mesh = inla.mesh.2d(loc= unique(data_est[, list(x,y,z)]),
                    cutoff=0.006,
                    min.angle=c(25,25),
                    max.edge=c(0.06,500) )
print(paste("New mesh constructed:", spatial_mesh$n, "vertices"))

# generate spde matern model from mesh
spde_matern =inla.spde2.matern(spatial_mesh,alpha=2) # formerly spde

# generate temporal mesh
temporal_mesh=inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2) # formerly mesh1d

# prep data for model fitting
cov_list<-data_est[, cov_names, with=F]
cov_list$year <- data_est$yearqtr

# TODO: work out whether gaul code and "empmean" covs shoudl be added here
cov_list <-as.list(cov_list)


# generate observation matrix
A_est =
  inla.spde.make.A(spatial_mesh, 
                   loc=as.matrix(data_est[, list(x,y,z)]), 
                   group=data_est$yearqtr,
                   group.mesh=temporal_mesh)
field_indices = inla.spde.make.index("field", n.spde=spatial_mesh$n,n.group=temporal_mesh$m)

# what is this?
stack_est = inla.stack(data=list(response=data_est$ihs_emp_access_dev),
                       A=list(A_est,1),
                       effects=
                         list(c(field_indices,
                                list(Intercept=1)),
                              c(cov_list)),
                       tag="est", remove.unused=TRUE)
stack_est<-inla.stack(stack_est)

model_formula<- as.formula(paste(
                          paste("response ~ -1 + Intercept  + "),
                          paste("f(field, model=spde_matern, group=field.group, control.group=list(model='ar1')) + ",sep=""),
                          paste(cov_names,collapse='+'),
                          sep=""))

#-- Call INLA and get results --#
mod_pred =   inla(model_formula,
                  data=inla.stack.data(stack_est),
                  family=c("gaussian"),
                  control.predictor=list(A=inla.stack.A(stack_est), compute=TRUE,quantiles=NULL),
                  control.compute=list(cpo=TRUE,waic=TRUE),
                  keep=FALSE, verbose=TRUE,
                  control.inla= list(strategy = 'gaussian',
                                     int.strategy='ccd', # close composite design ?
                                     verbose=TRUE,
                                     step.factor=1,
                                     stupid.search=FALSE)
)

print(paste("Saving outputs to", output_fname))
# save(mod_pred, file=output_fname)
save.image(output_fname)
# inla.ks.plot(mod.pred$cpo$pit, punif)




