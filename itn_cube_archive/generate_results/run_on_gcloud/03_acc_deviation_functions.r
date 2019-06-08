require('raster')
require('rgdal')
library(INLA)
library(RColorBrewer)
library(cvTools)
library(zoo)
library(boot)

convert.nodes<-function(r){
  v<-getValues(r)
  NAs<-is.na(v)
  v<-v[!NAs]
  un<-unique(v)
  l=length(un)
  tmp=v
  for(i in 1:l){
    tmp[v==un[i]]=i	
  }
  v=tmp
  r[!NAs]=v
  return(r)
}

ll.to.xyz<-function(ll){

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

temporalInfo<- function(f){
  library(zoo)
  start_dates<-as.yearmon(paste(f[, "year_start"], "-", f[, "month_start"],sep=""))
  end_dates<-as.yearmon(paste(f[, "year_end"], "-", f[, "month_end"],sep=""))
  mid_dates<-(start_dates+end_dates)/2
  return(as.numeric(mid_dates))
}
load.data<-function(loc){
  f<-read.csv(loc)
  year=temporalInfo(f)	# get months
  # for spatial only data snap points in 5km pixel to centroid and then aggregate number examined and positive
  r<-raster('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif')
  NAvalue(r)=-9999
  ce<-cellFromXY(r,f[,c('lon','lat')])  # get cell numbers
  cexy<-xyFromCell(r,ce)	# get xy from cells
  data<-as.data.frame(f)
  data$lon<-cexy[,1]
  data$lat<-cexy[,2]
  data$year<-floor(year)
  xyz<-ll.to.xyz(data[,c('lon','lat')])
  data$x<-xyz[,'x']
  data$y<-xyz[,'y']
  data$z<-xyz[,'z']
  return(data)
}

reposition.points<-function(ls,points,square){
  
  rm<-extract(ls,cbind(points[,'lon'],points[,'lat']))	#get all point values
  index<-which(is.na(rm))
  corindex<-c()
  
  for(x in index){
    row1 =  rowFromY(ls,points[x,'lat'])-square 
    col1 = colFromX(ls,points[x,'lon'])-square
    if(!is.na(col1) & !is.na(row1)){ #check for invalide boundaries
      mat<-raster::getValuesBlock(ls, row=row1, nrows=square*2+1, col=col1, ncols=square*2+1)	# get square around point
      mat <- matrix(mat, nrow=square*2+1, ncol=square*2+1, byrow=TRUE)
      matcells<-mat
      for(i in 1:dim(mat)[1]){
        for(j in 1:dim(mat)[2]){
          matcells[i,j]=as.numeric(cellFromRowCol(ls,(row1-1+i),(col1-1+j)))	# get cell numbers
        }
      }
      matcoords<-xyFromCell(ls,matcells) # get raw lat lons
      matcells<-as.vector(matcells)
      if(any(is.na(mat)==FALSE)){	#is any are not NA
        matcells<-matcells[!is.na(as.vector(mat))]
        matcoords<-matcoords[!is.na(as.vector(mat)),]
        dist<-pointDistance(cbind(points[x,'lon'],points[x,'lat']),matcoords,longlat=T)	# find the closest non na valid cell
        newp<-xyFromCell(ls,matcells[which.min(dist)])
        points[x,'lon']<-newp[,'x']
        points[x,'lat']<-newp[,'y']
        st<-paste('<',points$id[x],'-',points$name[x],'>')
        print(paste('-->Point',st ,'has been repositioned at distance of',min(dist)/1000,"km","at(",points[x,'lon'],",",points[x,'lat'],")"))  	
        corindex<-c(corindex,x)
      } else {
        st<-paste('<',points$id[x],'-',points$name[x],'>')	  
        print(paste("-->Point ",st,points[x,'lon'],points[x,'lat'],"is not positioned properly"))
      }
    } else {
      print(paste('Invalid boundry at',points[x,'lon'],points[x,'lat']))
    }
  }
  return(points)
  
}

ll2xyz <- function (longlat, radius = 6371) {
  
  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(longlat, 'matrix') | inherits(longlat, 'data.frame'))
  stopifnot(ncol(longlat) == 2)
  
  # extract required columns
  longitude <- longlat[, 1]
  latitude <- longlat[, 2]
  
  # convert
  ans <- data.frame(x = radius * cos(latitude) * cos(longitude),
                    y = radius * cos(latitude) * sin(longitude),
                    z = radius * sin(latitude))
  
  return (ans)
  
}



xyz2ll <- function (xyz, radius = 6371) {
  
  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(longlat, 'matrix') | inherits(longlat, 'data.frame'))
  stopifnot(ncol(xyz) == 3)
  
  # extract columns
  x <- xyz[, 1]
  y <- xyz[, 2]
  z <- xyz[, 3]
  
  # convert
  ans <- data.frame(longitude = atan2(y, x),
                    latitude = asin(z / radius))
  
  return (ans)
  
}



emplogit<-function(Y,N) {top=Y*N+0.5;bottom=N*(1-Y)+0.5;return(log(top/bottom))} # empirical logit number
emplogit2<-function(Y,N) {top=Y+0.5;bottom=N-Y+0.5;return(log(top/bottom))} # empirical logit two numbers


#from data load block:
IHS <- function(x, theta){  # Inverse IHS transformation
  (1/theta)*asinh(theta * x)
}

Inv.IHS <- function(x, theta){  # IHS transformation
  (1/theta)*sinh(theta * x)
}


IHS.loglik <- function(theta,x){
  
  IHS <- function(x, theta){  # function to IHS transform
    asinh(theta * x)/theta
  }
  
  n <- length(x)
  xt <- IHS(x, theta)
  
  log.lik <- -n*log(sum((xt - mean(xt))^2))- sum(log(1+theta^2*x^2))
  return(log.lik)
}


