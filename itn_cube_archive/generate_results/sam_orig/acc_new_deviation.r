require('raster')
require('rgdal')
library(INLA)
library(RColorBrewer)
library(cvTools)
library(zoo)
library(boot)
source('/home/backup/Space Time Malaria/algorithm.V1.R')
source('/home/backup/Space Time Malaria/INLAFunctions.R')

################# function block ######################################
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



########################################################################################################################
##### load point #######  ##########################################################################################

data<-read.csv('/home/backup/ITNcube/ITN_final_clean_access_20thNov2017.csv')
data<-data[!is.na(data$year),]

##########################################################################################
# covariate block
library(zoo)
foldersd<-c('/home/drive/cubes/5km/LST_day/mean/',
#'/home/drive/cubes/5km/LST_delta/mean/',
'/home/drive/cubes/5km/LST_night/mean/',
'/home/drive/cubes/5km/EVI/mean/',
#'/home/drive/cubes/5km/TCB/mean/',
'/home/drive/cubes/5km/TCW/mean/',
'/home/drive/cubes/5km/TSI/mean/')


dates<-as.Date(as.yearmon(data$year))
library(stringr)
split.dates<-str_split_fixed(dates, "-", 3)
uniques<-unique(split.dates)
l=length(foldersd)
cell<-data$cellnumber

#compute covariates for unique year-months		
library(doParallel)
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
#assign these to the correct cells in data
dynamic.covs<-matrix(NA,nrow=nrow(data),ncol=length(foldersd))
for(i in 1:nrow(uniques)){
	wh<-split.dates[,1]==uniques[i,1] & split.dates[,2]==uniques[i,2]
	dynamic.covs[wh,]<-covs.list.dyn[[i]]
}
#colnames(dynamic.covs)<-c('lst_day','lst_delta','lst_night','evy','tcb','tcw','tsi')
colnames(dynamic.covs)<-c('lst_day','lst_night','evy','tcw','tsi')

######## year only

foldery<-c('/home/drive/cubes/5km/IGBP_Landcover/Fraction/','/home/drive/cubes/5km/AfriPop/')

uniques<-unique(split.dates[,1])


l=17+1 # 17 fraction classes and 1 afripop
cell<-data$cellnumber

#compute covariates for unique year-months		
library(doParallel)
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

#### static covariates
#### folders
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



##########################################################################################
#data load block
require('raster')
require('rgdal')
library(INLA)
library(RColorBrewer)
library(cvTools)
library(zoo)
library(boot)

library(dismo)
library(gbm)
########## more functions
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
##############

load('/home/backup/ITNcube/preload.Rdata')

data<-data[complete.cases(all.covs),]
all.covs<-all.covs[complete.cases(all.covs),]

all.covs<-as.data.frame(all.covs)


covariate.names<-colnames(all.covs)


data<-cbind(data,all.covs)

data<-data[sample(1:nrow(data),replace=F),]




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

##########################################################################################
#indicator block
drow<-nrow(data) ## Number of data points
# load country rasters
cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep=""))
NAvalue(cn)<--9999
POPULATIONS<-read.csv('/home/backup/ITNcube/country_table_populations.csv') # load table to match gaul codes to country names

data$yearqtr[data$yearqtr>=2015]=2014.75


data$accdev<-emplogit2(data$P,data$N)-emplogit(data$Amean,1000)
theta<-optimise(IHS.loglik, lower=0.001, upper=50, x=data$accdev, maximum=TRUE) # 2.0407
theta<-theta$maximum
data$accdev=IHS(data$accdev,theta)

formula<- as.formula(paste(paste("accdev ~ "),paste(covariate.names,collapse='+'),sep=""))
l<-lm(formula,data=data)

#library(mgcv)
#formula<- as.formula(paste(paste("accdev ~ "),paste('s(',covariate.names,')',collapse='+'),sep=""))
#l<-gam(formula,data=data)

lpacc<-plogis(emplogit(data$Amean,1000)+l$fitted.values)
#plot(data$P/data$N,lpacc)
cor(data$P/data$N,lpacc)
mean(abs((data$P/data$N-lpacc)))


INLA:::inla.dynload.workaround() 

xyz<-ll.to.xyz(cbind(data$lon,data$lat))
data<-cbind(data,xyz)

data.est<-data
un<-paste(data.est$x,data.est$y,data.est$z,sep=':')
dup<-!duplicated(un)
#mesh = inla.mesh.2d(loc=cbind(data.est[dup,'x'],data.est[dup,'y'],data.est[dup,'z']),
#	cutoff=0.006,
#	min.angle=c(25,25),
#	max.edge=c(0.015,500) )

mesh = inla.mesh.2d(loc=cbind(data.est[dup,'x'],data.est[dup,'y'],data.est[dup,'z']),
	cutoff=0.006,
	min.angle=c(25,25),
	max.edge=c(0.06,500) )





print(paste('this is the number of mesh verticies',mesh$n))

spde=inla.spde2.matern(mesh,alpha=2)

mesh1d=inla.mesh.1d(seq(2000,2015,by=2),interval=c(2000,2015),degree=2)

data.cov<-data.est[,covariate.names]
est.cov<-as.list(data.cov)
est.cov$year<-data.est$yearqtr
est.cov$gaul<-data$gaul
est.cov$empmean<-data$empmean
A.est =
inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'yearqtr'],group.mesh=mesh1d)
field.indices = inla.spde.make.index("field", n.spde=mesh$n,n.group=mesh1d$m)

stack.est = inla.stack(data=list(response=data.est$accdev),
				A=list(A.est,1),
				effects=
				list(c(field.indices,
				list(Intercept=1)),
				c(est.cov)),
				tag="est", remove.unused=TRUE)

				


#
#formula1<- as.formula(paste(
#paste("response ~ -1  + Intercept  "),
#sep=""))
#
#formula2<- as.formula(paste(
#paste("response ~ -1 + Intercept  + "),
#paste(covariate.names,collapse='+'),
#sep=""))
#
#formula3<- as.formula(paste(
#paste("response ~ -1 + Intercept  + "),
#paste("f(field, model=spde) + ",sep=""),
#paste(covariate.names,collapse='+'),
#sep=""))
#
#formula4<- as.formula(paste(
#paste("response ~ -1 + Intercept  + "),
#paste("f(field, model=spde,group=field.group, control.group=list(model='ar1')) + ",sep=""),
#paste(covariate.names,collapse='+'),
#sep=""))
#
#
#formulas<-list(formula1,formula2,formula3,formula4)
#
#
#stack.est<-inla.stack(stack.est)
#library(doMC)
#registerDoMC(5)
#models<-foreach(i=1:4) %dopar% {
#
#	formula<-formulas[[i]]
#	#-- Call INLA and get results --#
#	mod.pred =   inla(formula,
#				 data=inla.stack.data(stack.est),
#				 family=c("gaussian"),
#				 ##################################################
#				 ##################################################
#				 control.predictor=list(A=inla.stack.A(stack.est), compute=TRUE,quantiles=NULL),
#				 control.compute=list(cpo=TRUE,waic=TRUE),
#				 keep=FALSE, verbose=TRUE,
#				 control.inla= list(strategy = 'gaussian',
#				 int.strategy='ccd',
#				 verbose=TRUE,
#				 step.factor=1,
#				 mode.known=FALSE,
#				 stupid.search=FALSE)
#		 )      
#		return(mod.pred)
#}

formula1<- as.formula(paste(
paste("response ~ -1 + Intercept  + "),
paste("f(field, model=spde,group=field.group, control.group=list(model='ar1')) + ",sep=""),
paste(covariate.names,collapse='+'),
sep=""))

stack.est<-inla.stack(stack.est)

#-- Call INLA and get results --#
mod.pred =   inla(formula1,
		 data=inla.stack.data(stack.est),
		 family=c("gaussian"),
		 control.predictor=list(A=inla.stack.A(stack.est), compute=TRUE,quantiles=NULL),
		 control.compute=list(cpo=TRUE,waic=TRUE),
		 keep=FALSE, verbose=TRUE,
		 control.inla= list(strategy = 'gaussian',
		 int.strategy='ccd',
		 verbose=TRUE,
		 step.factor=1,
		 stupid.search=FALSE)
 )      



save.image('/home/backup/ITNcube/ITN_cube_access_dynamic_access deviation_21112017.Rdata')

inla.ks.plot(mod.pred$cpo$pit, punif)

waic<-c()
for(i in 1:4){
 waic[i]<-models[[i]]$waic$waic
}
a=(waic[1]-waic[2])/waic[1]
b=(waic[2]-waic[3])/waic[1]
c=(waic[3]-waic[4])/waic[1]

a/(a+b+c)
b/(a+b+c)
c/(a+b+c)

#mod.pred=models[[4]]
index= inla.stack.index(stack.est,"est")$data
lp=mod.pred$summary.linear.predictor$mean[index]
lp=Inv.IHS(lp,theta)

lpgap<-plogis(emplogit(data$Amean,1000)+lp)
plot(data$P/data$N,lpgap)
cor(data$P/data$N,lpgap)
mean(abs((data$P/data$N-lpgap)))

a=mod.pred$summary.fixed$mean
sig<-sign(mod.pred$summary.fixed$'0.025quant')==sign(mod.pred$summary.fixed$'0.975quant')
names(a)<-rownames(mod.pred$summary.fixed)
sort(a[sig])

##### cross validation
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
                 
                                

## Space-time precision: the magic
Q = kronecker(Q.time, Q.space)

#A =inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']))
A= inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'yearqtr'],group.mesh=mesh1d)

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

#plot(mod.pred$summary.linear.predictor$sd[index],(sqrt(V)))

mval<-rep(NA,nrow(data))
for(i in 1:nrow(data)){
	mval[i]<-y[i]-(y[i]-mu_xy[i])/(1-(q[i]*V[i]))
}
mval<-Inv.IHS(mval,theta)
mval<-plogis(emplogit(data$Amean,1000)+mval)
plot(data$P/data$N,mval)
cor(data$P/data$N,mval)
mean(abs((data$P/data$N-mval)))




#annomally
plot(y,mval)
plot(data$P/data$N,plogis(emplogit(data$Amean,1000)+mval))
cor(data$P/data$N,plogis(emplogit(data$Amean,1000)+mval))
mean(abs(data$P/data$N-plogis(emplogit(data$Amean,1000)+mval)))


un=unique(data$Survey)
actual<-predicted<-c()

for(i in 1:length(un)){

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



#save.image('/home/backup/ITNcube/ITN_cube_access_dynamic.Rdata')
