require('raster')
require('rgdal')
library(INLA)
library(RColorBrewer)
library(cvTools)
library(zoo)
library(boot)
library(stringr)
library(doParallel)

source('/home/backup/Space Time Malaria/algorithm.V1.R')
source('/home/backup/Space Time Malaria/INLAFunctions.R')
source("acc_deviation_functions.r")

# Load data from create_database.r  ------------------------------------------------------------
data<-read.csv('/home/backup/ITNcube/ITN_final_clean_access_20thNov2017.csv')
data<-data[!is.na(data$year),]

# Load dynamic covariates  ------------------------------------------------------------
foldersd<-c('/home/drive/cubes/5km/LST_day/mean/',
            #'/home/drive/cubes/5km/LST_delta/mean/',
            '/home/drive/cubes/5km/LST_night/mean/',
            '/home/drive/cubes/5km/EVI/mean/',
            #'/home/drive/cubes/5km/TCB/mean/',
            '/home/drive/cubes/5km/TCW/mean/',
            '/home/drive/cubes/5km/TSI/mean/')

# find raster cell ids from access/use data 
dates<-as.Date(as.yearmon(data$year))
split.dates<-str_split_fixed(dates, "-", 3)
uniques<-unique(split.dates)
l=length(foldersd)
cell<-data$cellnumber

# extract covariates for unique year-months		
registerDoParallel(62)
covs.list.dyn<-foreach(i=1:nrow(uniques)) %dopar% { # loop through unique names
  wh<-split.dates[,1]==uniques[i,1] & split.dates[,2]==uniques[i,2]
  un.cells<-cell[wh]
  tmp<-matrix(NA,nrow=length(un.cells),ncol=l)
  for(j in 1:l){
    year<-as.numeric(uniques[i,1])
    
    while(!file.exists(paste0(foldersd[j],year,'.',uniques[i,2],'.mean.tif'))){	# some years have no data yet
      year=year-1
    }
    r=raster(paste0(foldersd[j],year,'.',uniques[i,2],'.mean.tif'))
    NAvalue(r)=-9999
    tmp[,j]<-r[un.cells]
  }
  return(tmp)
}

# assign these to the correct cells in data
dynamic.covs<-matrix(NA,nrow=nrow(data),ncol=length(foldersd))
for(i in 1:nrow(uniques)){
  wh<-split.dates[,1]==uniques[i,1] & split.dates[,2]==uniques[i,2]
  dynamic.covs[wh,]<-covs.list.dyn[[i]]
}
#colnames(dynamic.covs)<-c('lst_day','lst_delta','lst_night','evy','tcb','tcw','tsi')
colnames(dynamic.covs)<-c('lst_day','lst_night','evy','tcw','tsi')

# Load year-only (??) covariates  ------------------------------------------------------------
foldery<-c('/home/drive/cubes/5km/IGBP_Landcover/Fraction/','/home/drive/cubes/5km/AfriPop/')

uniques<-unique(split.dates[,1])

l=17+1 # 17 fraction classes and 1 afripop
cell<-data$cellnumber

#compute covariates for unique year-months	 (surely just years? --abv)	
# library(doParallel) # do you need to reload this library every time?
registerDoParallel(62)
covs.list.year<-foreach(i=1:length(uniques)) %dopar% { # loop through unique names
  wh<-split.dates[,1]==uniques[i] 
  un.cells<-cell[wh]
  tmp<-matrix(NA,nrow=length(un.cells),ncol=l)
  for(j in 1:(l-1)){
    # no land cover for 2013
    if(uniques[i]>2012){
      r=raster(paste0(foldery[1],'2012','.fraction.class.',j-1,'.tif'))
    }else{
      r=raster(paste0(foldery[1],uniques[i],'.fraction.class.',j-1,'.tif'))
    }		
    NAvalue(r)=-9999
    tmp[,j]<-r[un.cells]
  }
  if(uniques[i]>2015){
    r=raster(paste0(foldery[2],'2015','.total.population.tif'))
  }else{
    r=raster(paste0(foldery[2],uniques[i],'.total.population.tif'))
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


# Load static covariates  ------------------------------------------------------------

folderss<-c('/home/drive/cubes/5km/Topographic/Africa_TMI_90m.mean.tif',
            '/home/drive/cubes/5km/Topographic/Africa_SRTM_90m.mean.tif',
            '/home/drive/cubes/5km/Topographic/Africa_slope_90m.mean.tif',
            #'/home/drive/cubes/5km/Topographic/Africa_FA_90m.mean.tif',
            '/home/drive/cubes/5km/Seasonality/pf_seasonality.tif',
            '/home/drive/cubes/5km/Poverty/PET_1950-2000_5km.mean.tif',
            '/home/drive/cubes/5km/Poverty/AI_1950-2000_5km.mean.tif',
            '/home/drive/cubes/5km/Poverty/accessibility_50k_5km.mean.tif',
            #'/home/drive/cubes/5km/worldclim/prec57a0.tif',
            '/home/drive/cubes/5km/Poverty/viirs_nighttime_5km.mean.tif')
#'/home/drive/cubes/5km/Poverty/DMSP_F18_5km.mean.tif')


st<-stack(folderss)
NAvalue(st)=-9999

static.covs<-st[data$cellnumber]

all.covs<-cbind(static.covs,yearonly.covs,dynamic.covs)

save.image('/home/backup/ITNcube/preload.Rdata')




