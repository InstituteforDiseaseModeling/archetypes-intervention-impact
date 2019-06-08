###############################################################################################################
## 04_use_gap.r
## Amelia Bertozzi-Villa
## May 2019
## 
## A restructuring of Sam Bhatt"s original code to run the INLA model for use gap (difference between access
## and use.)
## 
##############################################################################################################

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam func_dir=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor --input CODE=gs://map_data_z/users/amelia/itn_cube/code/amelia_refactor/04_use_gap.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_replicate_sam/ --command 'Rscript ${CODE}'


rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam//"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190606_replicate_sam//"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir") 
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

# load relevant functions
source(file.path(func_dir, "03_05_general_functions.r"))

set.seed(212)

output_fname <- file.path(output_dir, "04_use_gap.Rdata")

## Load data 
data <- fread(file.path(input_dir, "02_data_covariates.csv"))
data <- data[order(row_id)]


## Prep for modeling ## ---------------------------------------------------------

# drop NAs
dropped_rows <- nrow(data) - nrow(data[complete.cases(data)])
print(paste("dropping", dropped_rows, "rows of data due to null values in covariates!"))
data <- data[complete.cases(data)]

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

## BUG?: (here and access dev: don't include country code in regression/don't merge s sudan and sudan)

## BUG?: here, year is capped at 2014 instead of 2015 (but temporal mesh goes through 2015)
# adjust timing of data? this feels important****
data[, yearqtr:=pmin(yearqtr, 2013.75)]


# calculate use gap for data points (USING GAP2 FROM DATABASE CREATION), 
# transform via empirical logit and inverse hyperbolic sine

theta<-optimise(ihs_loglik, lower=0.001, upper=50, x=data$gap2, maximum=TRUE) 
theta_use<-theta$maximum

data[, ihs_gap2:=ihs(gap2, theta_use)] 

# transform data from latlong to cartesian coordinates
xyz<-ll_to_xyz(data[, list(row_id, longitude=lon, latitude=lat)])

data <- merge(data, xyz, by="row_id", all=T)

# shuffle row order (why?)
data <- data[sample(1:nrow(data),replace=F),]

# initialize inla
INLA:::inla.dynload.workaround() 

# make a copy of the data for estimation
data_est <- copy(data)

# generate spatial mesh using unique xyz values 
# TODO: is this all xyz is used for?
spatial_mesh_use = inla.mesh.2d(loc= unique(data_est[, list(x,y,z)]),
                            cutoff=0.006,
                            min.angle=c(25,25),
                            max.edge=c(0.06,500) )
print(paste("New mesh constructed:", spatial_mesh_use$n, "vertices"))

# generate spde matern model from mesh
spde_matern =inla.spde2.matern(spatial_mesh_use,alpha=2) # formerly spde

# generate temporal mesh
temporal_mesh=inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2) # formerly mesh1d

# prep data for model fitting
cov_list<-data_est[, cov_names, with=F]
cov_list$year <- data_est$yearqtr

# TODO: work out whether gaul code and "empmean" covs shoudl be added here
cov_list <-as.list(cov_list)


# generate observation matrix
A_est =
  inla.spde.make.A(spatial_mesh_use, 
                   loc=as.matrix(data_est[, list(x,y,z)]), 
                   group=data_est$yearqtr,
                   group.mesh=temporal_mesh)
field_indices = inla.spde.make.index("field", n.spde=spatial_mesh_use$n,n.group=temporal_mesh$m)

# Generate "stack"
stack_est = inla.stack(data=list(response=data_est$ihs_gap2),
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
  paste(cov_names,collapse="+"),
  sep=""))

#-- Call INLA and get results --#
mod_pred_use =   inla(model_formula,
                  data=inla.stack.data(stack_est),
                  family=c("gaussian"),
                  control.predictor=list(A=inla.stack.A(stack_est), compute=TRUE,quantiles=NULL),
                  control.compute=list(cpo=TRUE,waic=TRUE),
                  keep=FALSE, verbose=TRUE,
                  control.inla= list(strategy = "gaussian",
                                     int.strategy="ccd", # close composite design ?
                                     verbose=TRUE,
                                     step.factor=1,
                                     stupid.search=FALSE)
)

print(summary(mod_pred_use))

print(paste("Saving outputs to", output_fname))
save(mod_pred_use, spatial_mesh_use, theta_use, file=output_fname)
# save.image(output_fname)




