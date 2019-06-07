
# STEP 2: COVARIATES

# unused covariates:
# '/home/drive/cubes/5km/LST_delta/mean/' #
# '/home/drive/cubes/5km/TCB/mean/' #
# '/home/drive/cubes/5km/Topographic/Africa_FA_90m.mean.tif' #
# '/home/drive/cubes/5km/worldclim/prec57a0.tif'


# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/cubes_5km func_dir=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud --input database_fname=gs://map_data_z/users/amelia/itn_cube/results/20190606_rerun_sam/01_database.csv CODE=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud/02_generate_covariates.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_rerun_sam/ --command 'Rscript ${CODE}'

rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "doParallel", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "stringr", "dismo", "gbm"))

if(Sys.getenv("input_dir")=="") {
  database_fname <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/output/ITN_final_clean_access_9Feb2019.csv"
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/input"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/access_deviation"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  database_fname <- Sys.getenv("database_fname") # location of output file from generate_database_refactored.r
  input_dir <- Sys.getenv("input_dir") # here, location of covariate data 
  output_dir <- Sys.getenv("output_dir") 
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

# load relevant functions
# source(file.path(func_dir, "acc_deviation_functions.r"))
# source(file.path(func_dir, "algorithm.V1.R"))
# source(file.path(func_dir, "INLAFunctions.R"))


output_fname <- file.path(output_dir, "02_covariates.Rdata")

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

save.image(output_fname)
