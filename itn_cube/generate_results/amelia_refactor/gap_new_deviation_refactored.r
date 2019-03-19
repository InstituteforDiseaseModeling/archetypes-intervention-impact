
##########################################################################################


# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-16 --logging gs://map_data_z/users/amelia/logs --input-recursive joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data func_dir=gs://map_data_z/users/amelia/itn_cube/code --input CODE=gs://map_data_z/users/amelia/itn_cube/code/gap_new_deviation_refactored.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/use_gap --command 'Rscript ${CODE}'


rm(list=ls())

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "data.table", "rgdal", "INLA", "RColorBrewer", "cvTools", "boot", "dismo", "gbm"))

if(Sys.getenv("joint_dir")=="") {
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/use_gap"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  output_dir <- Sys.getenv("output_dir")
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

# # load relevant functions
# source(file.path(func_dir, "acc_deviation_functions.r"))
# source(file.path(func_dir, "algorithm.V1.R"))
# source(file.path(func_dir, "INLAFunctions.R"))
# 
# 
load(file.path(joint_dir, 'preload.Rdata')) # same as access_new_deviation.R


## temp while you still load the entire environment with 'preload' (todo: save only the necessary components of the image) --------
if(Sys.getenv("joint_dir")=="") {
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/use_gap"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  joint_dir <- Sys.getenv("joint_dir") # location for shared datasets across itn cube scripts
  output_dir <- Sys.getenv("output_dir")
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}
##------


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
cn<-raster(file.path(joint_dir, 'african_cn5km_2013_no_disputes.tif'))
NAvalue(cn)<--9999
POPULATIONS<-read.csv(file.path(joint_dir, 'country_table_populations.csv')) # load table to match gaul codes to country names

#### change south sudan to sudan
#POPULATIONS[POPULATIONS$NAME=='South Sudan','NAME']='Sudan'
gauls<-extract(cn,data[,c('lon','lat')]) # extract gaul code for points
data$african_cn<-gauls #update the gaul codes
rm(gauls)


data$yearqtr[data$yearqtr>=2014]=2013.75

data<-data[!is.na(data$gap2),]

theta<-optimise(IHS.loglik, lower=0.001, upper=50, x=data$gap2, maximum=TRUE) # 2.0407
theta<-theta$maximum
data$gapIHS=IHS(data$gap2,theta)

formula<- as.formula(paste(paste("gapIHS ~ "),paste(covariate.names,collapse='+'),sep=""))
l<-lm(formula,data=data)

plot(l$fitted.values,data$gapIHS)
cor(l$fitted.values,data$gapIHS)
abline(0,1,col='red')


lpgap<-plogis(emplogit2(data$P,data$N)-Inv.IHS(l$fitted.values,theta))
plot(data$Pu/data$N,lpgap)
cor(data$Pu/data$N,lpgap)
mean(abs((data$Pu/data$N-lpgap)))



xyz<-ll.to.xyz(cbind(data$lon,data$lat))
data<-cbind(data,xyz)


INLA:::inla.dynload.workaround() 

data.est<-data
un<-paste(data.est$x,data.est$y,data.est$z,sep=':')
dup<-!duplicated(un)
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

stack.est = inla.stack(data=list(response=data.est$gapIHS),
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


# save.image(file.path(output_dir, 'ITN_cube_gap_dynamic_18March2019.Rdata.Rdata'))

# inla.ks.plot(models[[4]]$cpo$pit, punif)
# 
# waic<-c()
# for(i in 1:4){
#  waic[i]<-models[[i]]$waic$waic
# }
# a=(waic[1]-waic[2])/waic[1]
# b=(waic[2]-waic[3])/waic[1]
# c=(waic[3]-waic[4])/waic[1]
# 
# a/(a+b+c)
# b/(a+b+c)
# c/(a+b+c)
# 
# mod.pred=models[[4]]
index= inla.stack.index(stack.est,"est")$data
lp=mod.pred$summary.linear.predictor$mean[index]
lp=Inv.IHS(lp,theta)

lpgap<-plogis(emplogit2(data$P,data$N)-lp)
plot(data$Pu/data$N,lpgap)
cor(data$Pu/data$N,lpgap)
mean(abs((data$Pu/data$N-lpgap)))

a=mod.pred$summary.fixed$mean
sig<-sign(mod.pred$summary.fixed$'0.025quant')==sign(mod.pred$summary.fixed$'0.975quant')
names(a)<-rownames(mod.pred$summary.fixed)
sort(a[sig])

##### cross validation
y=data.est$gapIHS

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
A= inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'flooryear'],group.mesh=mesh1d)

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

save.image(file.path(output_dir, 'ITN_cube_gap_dynamic_18March2019.Rdata.Rdata'))
