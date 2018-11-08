
library(VGAM)
#set window
library(maptools)

library(rgeos)

library(png)

library(sp)

World<-readShapePoly('/home/backup/MAP_templates/Africa_area/afri_word.shp')

World <- gSimplify(World, tol=0.1, topologyPreserve=TRUE)

Africa<-readShapePoly('/home/backup/MAP_templates/Africa/Africa.shp')

Africa <- gSimplify(Africa, tol=0.1, topologyPreserve=TRUE)

Water<-readShapePoly('/home/backup/MAP_templates/Africa_area/africa_water.shp')

Water <- gSimplify(Water, tol=0.1, topologyPreserve=TRUE)

img <- readPNG('/home/backup/MAP_templates/Logo/map-logo_bg.png')



require('raster')
require('rgdal')
library(INLA)
library(RColorBrewer)
library(zoo) 
#load('/home/backup/ITNcube/ITN_cube_access_dynamic_access deviation.Rdata')
#load('/home/backup/ITNcube/ITN_cube_access_dynamic_access deviation_nocovariatetransform.Rdata')
load('/home/backup/ITNcube/ITN_cube_access_dynamic_access deviation_21112017.Rdata')
#mod.pred.acc=models[[4]]
mod.pred.acc=mod.pred

theta.acc<-theta
mesh.acc<-mesh

load('/home/backup/ITNcube/ITN_cube_gap_dynamic_21112017.Rdata')
mod.pred.use=mod.pred
theta.use<-theta
mesh.use<-mesh


#Create Map template
create_access_cover_mat<-function(){
	nr=40
	access_cover_mat<-matrix(data=FALSE,nrow=nr,ncol=10)
	index=1:ncol(access_cover_mat)
	for(i in 1:nrow(access_cover_mat)){
		access_cover_mat[i,index<=i*2]=TRUE
	}
	access_cover_mat<-rbind(rep(FALSE,10),access_cover_mat)

	return(access_cover_mat)
	}

calc_access<-function(accdate1,func0,func1,hh,nc=10,nr=40){

		k1<-k0<-rep(NA,10)
		for(j in 1:length(func0)){
			k1[j]<-func1[[j]](accdate1)
			k0[j]<-func0[[j]](accdate1)
		}
		k0[is.nan(k0)]<-0
		matModel<-matrix(data=0,nrow=nr+1,ncol=nc)
		cnames<-colnames(matModel)<-1:nc
		rnames<-rownames(matModel)<-0:nr
		
		matModel[1,]<-hh*k0

		remaining<-hh*(1-k0)
		for(i in 1:ncol(matModel)){
			l<-k1[i]
			matModel[2:nrow(matModel),i]<-dpospois(1:(nrow(matModel)-1),l)*remaining[i]
		}	
		access_cover_mat<-create_access_cover_mat()

		allind<-get.indicators.model(matModel,nc,nr)	
		return(allind[3])
}


get.indicators.model<-function(mat,nc,nr){
	cnames<-colnames(mat)<-1:nc
	rnames<-rownames(mat)<-0:nr	
	# formulas from kilians paper
	####### INDICATOR 1 ##########
	ind1<-sum(mat[rownames(mat)>0,])

	####### INDICATOR 2 ##########
	ind2mat<-matrix(data=FALSE,nrow=nrow(mat),ncol=ncol(mat))
	#detemrmine columns that meet 1 net between two people
	for(x in 1:ncol(mat)){
		 		ind2mat[,x]<- (2*rnames/x)>=1
	}	
	ind2<- sum(ind2mat*mat)

	####### INDICATOR 3 ##########
	ind3mat<-matrix(data=0,nrow=nrow(mat),ncol=ncol(mat))
	ind3tot<-matrix(data=0,nrow=nrow(mat),ncol=ncol(mat))

	for(x in 1:ncol(mat)){
		for(y in 1:nrow(mat)){
			ind3mat[y,x]<-2*(y-1)*mat[y,x] 
			ind3tot[y,x]<-x*mat[y,x] # total people
			if(ind3mat[y,x]>ind3tot[y,x]) ind3mat[y,x]=ind3tot[y,x] #cap on estimate
			
		}
	}
	ind3<-sum(ind3mat)/sum(ind3tot)

	####### INDICATOR 4 ##########
	ind4<-1-(ind2/ind1)
	return(c(ind1,ind2,ind3,ind4))
}


get.pred.locs<-function(){
	cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep="")) #load raster
	NAvalue(cn)=-9999
	pred_val<-getValues(cn)#get values again
	w<-is.na(pred_val) #find NAs again
	index<-1:length(w) 
	index<-index[!w]
	pred_locs<-xyFromCell(cn,1:ncell(cn))  #get prediction locations
	pred_locs<-pred_locs[!w,] #remove NA cells
	############### Plotting anomoly
	colnames(pred_locs)<-c('longitude','latitude')
	pred_locs<-ll.to.xyz(pred_locs) #get locations
	return(pred_locs)
}

######## year only



get.pred.covariates.dynamic<-function(year,month){
	library(zoo)
	foldersd<-c('/home/drive/cubes/5km/LST_day/mean/',
	#'/home/drive/cubes/5km/LST_delta/mean/',
	'/home/drive/cubes/5km/LST_night/mean/',
	'/home/drive/cubes/5km/EVI/mean/',
	#'/home/drive/cubes/5km/TCB/mean/',
	'/home/drive/cubes/5km/TCW/mean/',
	'/home/drive/cubes/5km/TSI/mean/')

	l=length(foldersd)

	cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep="")) #load raster
	NAvalue(cn)=-9999
	pred_val<-getValues(cn)#get values again
	w<-is.na(pred_val) #find NAs again
	index<-1:length(w) 
	index<-index[!w]
	pred_locs<-xyFromCell(cn,1:ncell(cn))  #get prediction locations
	pred_locs<-pred_locs[!w,] #remove NA cells
	############### Plotting anomoly
	colnames(pred_locs)<-c('longitude','latitude')

	#compute covariates for unique year-months		
	library(doParallel)
	registerDoParallel(62)
	covs.list.dyn<-foreach(i=1:l,.combine='cbind') %dopar% { # loop through unique names
			r=raster(paste0(foldersd[i],year,'.',month,'.mean.tif'))
			NAvalue(r)=-9999
			return(r[index])
		}
	colnames(covs.list.dyn)<-c('lst_day','lst_night','evy','tcw','tsi')

	return(covs.list.dyn)
}

get.pred.covariates.yearonly<-function(year){
	foldery<-c('/home/drive/cubes/5km/IGBP_Landcover/Fraction/','/home/drive/cubes/5km/AfriPop/')
	cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep="")) #load raster
	NAvalue(cn)=-9999
	pred_val<-getValues(cn)#get values again
	w<-is.na(pred_val) #find NAs again
	index<-1:length(w) 
	index<-index[!w]
	pred_locs<-xyFromCell(cn,1:ncell(cn))  #get prediction locations
	pred_locs<-pred_locs[!w,] #remove NA cells
	############### Plotting anomoly
	colnames(pred_locs)<-c('longitude','latitude')
	l=17 # 17 fraction classes and 1 afripop

	library(doParallel)
	registerDoParallel(62)
	covs.list.year<-foreach(i=1:l,.combine='cbind') %dopar% { # loop through unique names
			# no land cover for 2013
			if(year>2012){
				r=raster(paste0(foldery[1],'2012','.fraction.class.',i-1,'.tif'))
			}else{
				r=raster(paste0(foldery[1],year,'.fraction.class.',i-1,'.tif'))
			}		
			NAvalue(r)=-9999
			return(r[index])
		}
	covs.list.year<-covs.list.year[,-14]	#  # remove landcover 13  - Urban and built-up   - for collinearity with population
	covnames<-0:16
	covnames<-covnames[-14]	
	r=raster(paste0(foldery[2],year,'.total.population.tif'))
	NAvalue(r)=-9999
	pop<-r[index]
	covs.list.year<-cbind(covs.list.year,pop)
	colnames(covs.list.year)<-c(paste0('landcover',covnames),'populatopn')

	return(covs.list.year)	
}

get.pred.covariates.static<-function(){

#### static covariates
	#### folders
	folderss<-c('/home/drive/cubes/5km/Topographic/Africa_TMI_90m.mean.tif',
	'/home/drive/cubes/5km/Topographic/Africa_SRTM_90m.mean.tif',
	'/home/drive/cubes/5km/Topographic/Africa_slope_90m.mean.tif',
#	'/home/drive/cubes/5km/Topographic/Africa_FA_90m.mean.tif',
	'/home/drive/cubes/5km/Seasonality/pf_seasonality.tif',
	'/home/drive/cubes/5km/Poverty/PET_1950-2000_5km.mean.tif',
	'/home/drive/cubes/5km/Poverty/AI_1950-2000_5km.mean.tif',
	'/home/drive/cubes/5km/Poverty/accessibility_50k_5km.mean.tif',
	#'/home/drive/cubes/5km/worldclim/prec57a0.tif',
	'/home/drive/cubes/5km/Poverty/viirs_nighttime_5km.mean.tif')
	#'/home/drive/cubes/5km/Poverty/DMSP_F18_5km.mean.tif')
	
	
	foldery<-c('/home/drive/cubes/5km/IGBP_Landcover/Fraction/','/home/drive/cubes/5km/AfriPop/')
	cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep="")) #load raster
	NAvalue(cn)=-9999
	pred_val<-getValues(cn)#get values again
	w<-is.na(pred_val) #find NAs again
	index<-1:length(w) 
	index<-index[!w]
	pred_locs<-xyFromCell(cn,1:ncell(cn))  #get prediction locations
	pred_locs<-pred_locs[!w,] #remove NA cells

	st<-stack(folderss)
	NAvalue(st)=-9999

	static.covs<-st[index]
}

time_points<-2000:2016

static.covs<-get.pred.covariates.static()


for(xxx in 1:length(time_points)){
	dics<-cors<-c()
	pred_year=time_points[xxx]
	pred_month1="01"
	pred_month2="02"
	pred_month3="03"
	pred_month4="04"
	pred_month5="05"
	pred_month6="06"
	pred_month7="07"
	pred_month8="08"
	pred_month9="09"
	pred_month10="10"
	pred_month11="11"
	pred_month12="12"
	pred_year_cov=pred_year
	if(time_points[xxx]<=2000) pred_year_cov=2001 # covariates are not full in 2000
	if(time_points[xxx]>=2013) pred_year_cov=2013 # covariates are not full in 2000

	pred_mat_y<-get.pred.covariates.yearonly(pred_year_cov)

####################################################################################################################
#extract covariates
	pred_mat1_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month1)
	pred_mat2_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month2)
	pred_mat3_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month3)
	pred_mat4_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month4)
	pred_mat5_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month5)
	pred_mat6_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month6)
	pred_mat7_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month7)
	pred_mat8_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month8)
	pred_mat9_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month9)
	pred_mat10_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month10)
	pred_mat11_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month11)
	pred_mat12_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month12)


	pred_mat1<-cbind(static.covs,pred_mat_y,pred_mat1_s)
	pred_mat2<-cbind(static.covs,pred_mat_y,pred_mat2_s)
	pred_mat3<-cbind(static.covs,pred_mat_y,pred_mat3_s)
	pred_mat4<-cbind(static.covs,pred_mat_y,pred_mat4_s)
	pred_mat5<-cbind(static.covs,pred_mat_y,pred_mat5_s)
	pred_mat6<-cbind(static.covs,pred_mat_y,pred_mat6_s)
	pred_mat7<-cbind(static.covs,pred_mat_y,pred_mat7_s)
	pred_mat8<-cbind(static.covs,pred_mat_y,pred_mat8_s)
	pred_mat9<-cbind(static.covs,pred_mat_y,pred_mat9_s)
	pred_mat10<-cbind(static.covs,pred_mat_y,pred_mat10_s)
	pred_mat11<-cbind(static.covs,pred_mat_y,pred_mat11_s)
	pred_mat12<-cbind(static.covs,pred_mat_y,pred_mat12_s)

####################################################################################################################
#standardise and transform covariates
	#for(i in 1:12){	
	#	for(j in 1:length(trans)) {
	#		eval(parse(text=paste0("pred_mat",i,"[,",j,"]<-t[[trans[",j,"]]](pred_mat",i,"[,",j,"])")))
	#		eval(parse(text=paste0("pred_mat",i,"[,",j,"]<-(pred_mat",i,"[,",j,"]-ms[",j,"])/sds[",j,"]")))
	#	}
	#}
####################################################################################################################
#get pred locations
	pred_locs<-get.pred.locs()
####################################################################################################################
#get field
    #########################################
	pred_year=time_points[xxx]	     ########
	pred_year_field=pred_year
	if(pred_year>=2015) pred_year_field=2015 # adjust prediction year for field
	
    #########################################
	INLA:::inla.dynload.workaround() 
	
	pred_year1=pred_year
	Apreds<-inla.spde.make.A(mesh.acc, loc=cbind(pred_locs[,'x'],pred_locs[,'y'],pred_locs[,'z']),group=rep(pred_year_field,nrow(pred_locs)),group.mesh=mesh1d)
	
	A.pred_dynamic1=Apreds



####################################################################################################################
	#### get access mean from stock and flow
	cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep="")) #load raster
	NAvalue(cn)=-9999
	pred_val<-getValues(cn)#get values again
	w<-is.na(pred_val) #find NAs again
	index<-1:length(w) 
	index<-index[!w]
	gauls<-cn[index]


	load('/home/backup/ITNcube/prop0prop1.rData')
	POPULATIONS<-read.csv('/home/backup/National_Config_Data.csv')
	Country.list<-Cout
	Country.gaul<-Cout
	for(i in 1:length(Country.gaul)){
		Country.gaul[i]<-POPULATIONS[as.character(POPULATIONS[,'ISO3'])==as.character(Country.list[i]),'GAUL_Code']
	}	
	times_ind<-seq(0.0,17,0.25)+2000
	Country.gaul = as.numeric(as.vector(Country.gaul))
	accdate1<-as.numeric(as.yearmon(paste(time_points[xxx],1,sep='-')))
	accdate2<-as.numeric(as.yearmon(paste(time_points[xxx],2,sep='-')))
	accdate3<-as.numeric(as.yearmon(paste(time_points[xxx],3,sep='-')))
	accdate4<-as.numeric(as.yearmon(paste(time_points[xxx],4,sep='-')))
	accdate5<-as.numeric(as.yearmon(paste(time_points[xxx],5,sep='-')))
	accdate6<-as.numeric(as.yearmon(paste(time_points[xxx],6,sep='-')))
	accdate7<-as.numeric(as.yearmon(paste(time_points[xxx],7,sep='-')))
	accdate8<-as.numeric(as.yearmon(paste(time_points[xxx],8,sep='-')))
	accdate9<-as.numeric(as.yearmon(paste(time_points[xxx],9,sep='-')))
	accdate10<-as.numeric(as.yearmon(paste(time_points[xxx],10,sep='-')))
	accdate11<-as.numeric(as.yearmon(paste(time_points[xxx],11,sep='-')))
	accdate12<-as.numeric(as.yearmon(paste(time_points[xxx],12,sep='-')))

	
	acc1_p<-rep(NA,length(gauls))
	acc2_p<-rep(NA,length(gauls))
	acc3_p<-rep(NA,length(gauls))
	acc4_p<-rep(NA,length(gauls))
	acc5_p<-rep(NA,length(gauls))
	acc6_p<-rep(NA,length(gauls))
	acc7_p<-rep(NA,length(gauls))
	acc8_p<-rep(NA,length(gauls))
	acc9_p<-rep(NA,length(gauls))
	acc10_p<-rep(NA,length(gauls))
	acc11_p<-rep(NA,length(gauls))
	acc12_p<-rep(NA,length(gauls))	

	for(i in 1:length(Country.list)){
		p<-out[[i]]
		p0<-p[,,1] # get p0 parameters
		p1<-p[,,2] # get p1 parameters	
		func1<-func0<-list()
		# make these into continuous functions
		for(f in 1:ncol(p0)){
			eval(parse(text=paste0('func0[[',f,']]<-approxfun(times_ind,p0[,',f,'])')))
		}
		for(f in 1:ncol(p0)){
			eval(parse(text=paste0('func1[[',f,']]<-approxfun(times_ind,p1[,',f,'])')))
		}
		# load household size distribution
		hh<-read.csv('/home/backup/Bucket model/HHsize.csv')
		KEY=read.csv('/home/backup/Bucket model/KEY.csv')
		cn_nm<-as.character(KEY[KEY$Name%in%Country.list[i],'Svy.Name'])
		hh_val<-hh[hh$HHSurvey%in%cn_nm,]
		if(nrow(hh_val)!=0){
			hh_val<-colSums(hh_val[,3:ncol(hh_val)])/sum(colSums(hh_val[,3:ncol(hh_val)]))
		} else {
			hh_val<-colSums(hh[,3:ncol(hh)])/sum(colSums(hh[,3:ncol(hh)]))
		}
		hh_size<-hh_val[1:10]
		hh_size[10]<-sum(hh_val[10:length(hh_val)])
		hh<-hh_size
		names(hh)<-1:10
		
		
	
		#### calculate access
	

		acc1<-calc_access(accdate1,func0,func1,hh)
		acc2<-calc_access(accdate2,func0,func1,hh)
		acc3<-calc_access(accdate3,func0,func1,hh)
		acc4<-calc_access(accdate4,func0,func1,hh)
		acc5<-calc_access(accdate5,func0,func1,hh)
		acc6<-calc_access(accdate6,func0,func1,hh)
		acc7<-calc_access(accdate7,func0,func1,hh)
		acc8<-calc_access(accdate8,func0,func1,hh)
		acc9<-calc_access(accdate9,func0,func1,hh)
		acc10<-calc_access(accdate10,func0,func1,hh)
		acc11<-calc_access(accdate11,func0,func1,hh)
		acc12<-calc_access(accdate12,func0,func1,hh)
		
		cc<-gauls==Country.gaul[i]
		acc1_p[cc]=emplogit(acc1,1000)	
		acc2_p[cc]=emplogit(acc2,1000)
		acc3_p[cc]=emplogit(acc3,1000)
		acc4_p[cc]=emplogit(acc4,1000)
		acc5_p[cc]=emplogit(acc5,1000)
		acc6_p[cc]=emplogit(acc6,1000)
		acc7_p[cc]=emplogit(acc7,1000)
		acc8_p[cc]=emplogit(acc8,1000)
		acc9_p[cc]=emplogit(acc9,1000)
		acc10_p[cc]=emplogit(acc10,1000)
		acc11_p[cc]=emplogit(acc11,1000)
		acc12_p[cc]=emplogit(acc12,1000)

	}
	##### put it all together
	
	Beta<-c(mod.pred.acc$summary.fixed[,1])
	Beta_cov<-Beta[2:length(Beta)]


	
	fixed_effects1<-Beta[1] + as.numeric(Beta_cov%*%t(pred_mat1))	
	fixed_effects2<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat2))	
	fixed_effects3<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat3))	
	fixed_effects4<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat4))	
	fixed_effects5<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat5))	
	fixed_effects6<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat6))	
	fixed_effects7<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat7))	
	fixed_effects8<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat8))	
	fixed_effects9<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat9))	
	fixed_effects10<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat10))	
	fixed_effects11<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat11))	
	fixed_effects12<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat12))	
	
			
	lp1=fixed_effects1 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp2=fixed_effects2 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp3=fixed_effects3 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp4=fixed_effects4 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp5=fixed_effects5 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp6=fixed_effects6 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp7=fixed_effects7 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp8=fixed_effects8 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp9=fixed_effects9 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp10=fixed_effects10 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp11=fixed_effects11 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
	lp12=fixed_effects12 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 


	lp1<-Inv.IHS(lp1,theta.acc)
	lp2<-Inv.IHS(lp2,theta.acc)
	lp3<-Inv.IHS(lp3,theta.acc)
	lp4<-Inv.IHS(lp4,theta.acc)
	lp5<-Inv.IHS(lp5,theta.acc)
	lp6<-Inv.IHS(lp6,theta.acc)
	lp7<-Inv.IHS(lp7,theta.acc)
	lp8<-Inv.IHS(lp8,theta.acc)
	lp9<-Inv.IHS(lp9,theta.acc)
	lp10<-Inv.IHS(lp10,theta.acc)
	lp11<-Inv.IHS(lp11,theta.acc)
	lp12<-Inv.IHS(lp12,theta.acc)


	lp1_lp = acc1_p+lp1
	lp2_lp = acc2_p+lp2
	lp3_lp = acc3_p+lp3
	lp4_lp = acc4_p+lp4
	lp5_lp = acc5_p+lp5
	lp6_lp = acc6_p+lp6
	lp7_lp = acc7_p+lp7
	lp8_lp = acc8_p+lp8
	lp9_lp = acc9_p+lp9
	lp10_lp = acc10_p+lp10
	lp11_lp = acc11_p+lp11
	lp12_lp = acc12_p+lp12
	
	acc_mean<-(acc1_p + acc2_p + acc3_p +acc4_p+acc5_p + acc6_p + acc7_p + acc8_p + acc9_p + acc10_p + acc11_p +acc12_p)/12
	
	lp = cbind(lp1,lp2,lp3,lp4,lp5,lp6,lp7,lp8,lp9,lp10,lp11,lp12) ## access deviation

	lp_lp = cbind(lp1_lp,lp2_lp,lp3_lp,lp4_lp,lp5_lp,lp6_lp,lp7_lp,lp8_lp,lp9_lp,lp10_lp,lp11_lp,lp12_lp) # access


	lp<-(lp)
	lp<-rowMeans(lp)

	lp_lp<-plogis(lp_lp)
	lp_lp<-rowMeans(lp_lp)


	colfunc <- colorRampPalette(c("blue","cyan","yellow","orange","red"))
	colfunc2 <- colorRampPalette(c("blue","cyan","yellow","orange","red"),bias=2)

	P<-cn
	P[!w]<-(lp_lp)
	P[!is.na(cn) & is.na(P)]=0

	P2<-cn
	P2[!w]<-(lp)
	P2[!is.na(cn) & is.na(P2)]=0

	P3<-cn
	P3[!w]<-plogis(acc_mean)
	P3[!is.na(cn) & is.na(P3)]=0

	writeRaster(P,paste0('/home/backup/ITNcube/ITN_',time_points[xxx],'.ACC.tif'),NAflag=-9999,overwrite=TRUE)
	writeRaster(P2,paste0('/home/backup/ITNcube/ITN_',time_points[xxx],'.DEV.tif'),NAflag=-9999,overwrite=TRUE)
	writeRaster(P3,paste0('/home/backup/ITNcube/ITN_',time_points[xxx],'.MEAN.tif'),NAflag=-9999,overwrite=TRUE)

#############################################################################################################################################
	#now predict use

	pred_year1=pred_year
	Apreds<-inla.spde.make.A(mesh.use, loc=cbind(pred_locs[,'x'],pred_locs[,'y'],pred_locs[,'z']),group=rep(pred_year,nrow(pred_locs)),group.mesh=mesh1d)
	
	A.pred_dynamic1=Apreds


	Beta<-c(mod.pred.use$summary.fixed[,1])
	Beta_cov<-Beta[2:length(Beta)]
	
			
	fixed_effects1<-Beta[1] + as.numeric(Beta_cov%*%t(pred_mat1))	
	fixed_effects2<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat2))	
	fixed_effects3<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat3))	
	fixed_effects4<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat4))	
	fixed_effects5<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat5))	
	fixed_effects6<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat6))	
	fixed_effects7<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat7))	
	fixed_effects8<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat8))	
	fixed_effects9<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat9))	
	fixed_effects10<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat10))	
	fixed_effects11<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat11))	
	fixed_effects12<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat12))	
	
			
	gap1=fixed_effects1 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap2=fixed_effects2 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap3=fixed_effects3 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap4=fixed_effects4 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap5=fixed_effects5 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap6=fixed_effects6 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap7=fixed_effects7 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap8=fixed_effects8 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap9=fixed_effects9 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap10=fixed_effects10 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap11=fixed_effects11 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	gap12=fixed_effects12 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
	

	gap1<-Inv.IHS(gap1,theta.use)
	gap2<-Inv.IHS(gap2,theta.use)
	gap3<-Inv.IHS(gap3,theta.use)
	gap4<-Inv.IHS(gap4,theta.use)
	gap5<-Inv.IHS(gap5,theta.use)
	gap6<-Inv.IHS(gap6,theta.use)
	gap7<-Inv.IHS(gap7,theta.use)
	gap8<-Inv.IHS(gap8,theta.use)
	gap9<-Inv.IHS(gap9,theta.use)
	gap10<-Inv.IHS(gap10,theta.use)
	gap11<-Inv.IHS(gap11,theta.use)
	gap12<-Inv.IHS(gap12,theta.use)

	use1<-(lp1_lp-gap1)
	use2<-(lp2_lp-gap2)
	use3<-(lp3_lp-gap3)
	use4<-(lp4_lp-gap4)
	use5<-(lp5_lp-gap5)
	use6<-(lp6_lp-gap6)
	use7<-(lp7_lp-gap7)
	use8<-(lp8_lp-gap8)
	use9<-(lp9_lp-gap9)
	use10<-(lp10_lp-gap10)
	use11<-(lp11_lp-gap11)
	use12<-(lp12_lp-gap12)
	
	lp_gap = cbind(gap1,gap2,gap3,gap4,gap5,gap6,gap7,gap8,gap9,gap10,gap11,gap12) # access
	lp_gap<-rowMeans(lp_gap)

	
	lp_use = cbind(use1,use2,use3,use4,use5,use6,use7,use8,use9,use10,use11,use12) # access
	lp_use<-plogis(lp_use)
	lp_use<-rowMeans(lp_use)
	# predict
	P4<-cn
	P4[!w]<-(lp_gap)
	P4[!is.na(cn) & is.na(P4)]=0

	P5<-cn
	P5[!w]<-(lp_use)
	P5[!is.na(cn) & is.na(P5)]=0

	writeRaster(P4,paste0('/home/backup/ITNcube/ITN_',time_points[xxx],'.GAP.tif'),NAflag=-9999,overwrite=TRUE)
	writeRaster(P5,paste0('/home/backup/ITNcube/ITN_',time_points[xxx],'.USE.tif'),NAflag=-9999,overwrite=TRUE)


}
##########################################################################################

library(raster)

library(VGAM)
#set window
library(maptools)

library(rgeos)

library(png)

library(sp)
time_points<-2000:2016
colfunc <- colorRampPalette(c("red","orange","yellow",'green'))

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))




match.cols<-function(val,n,colfunc){
	col<-data.frame(val=seq(min(val),max(val),length.out=n),col=colfunc(n))
	out<-rep(NA,length(col))
	for(i in 1:length(val)){
		out[i]<-as.character(col[which.min(abs(col$val-val[i])),'col'])
	}	
	return(out)
}


Plot_raster<-function(r,bias,title){
 n = 1000
 jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),bias=1.5) 
  lon_west=-17
  lon_east=54
  lat_south=-36.6
  lat_north=38
  par(bg="white")
  
  
  plot(1,1,type='n',xlim=c(lon_west,lon_east),ylim=c(lat_south,lat_north),xlab='',ylab='',axes=F,main=title)

  u <- par("usr")
  
  rect(u[1], u[3], u[2], u[4], col = rgb(180,225,255,maxColorValue=255), border = "black") 
  plot(World,add=T,col='white')

#### my raster
  image(r,add=T,legend=F,col=jet.colors(n))
  plot(Africa,add=T)
  plot(Water,col = rgb(180,225,255,maxColorValue=255),add=T) 
  r.range <- c(minValue(r), maxValue(r))
  leg.locs<-signif(seq(r.range[1],r.range[2],length.out=100),2)
  cols<-match.cols(leg.locs,n,jet.colors)
  label<-c(as.character(signif(r.range[1],2)),as.character(signif(r.range[2],2)))
  library(SDMTools)
  # xmin ymin xmax ymax
  rect(-17, -13, -8.5, -37, col = "white", border = "black") 
  pnts = cbind(x =c(-16.5,-12.5,-12.5,-16.5), y =c(-15,-15,-35.6,-35.6))
  legend.gradient(pnts,cols,label,title="")
  box()
}

Plot_points<-function(data,year,wh,title){
 data<-data[complete.cases(data),]
 if(wh=='use'){
 	val<-data$Pu/data$N
 } else {
 	val<-data$P/data$N
 
 }
 tmp_data<-data.frame(X=data$lon,Y=data$lat,T=data$flooryear,P=val)
 tmp_data<-rbind(tmp_data,c(-180,-180,year,1),c(-180,-180,year,0))
 tmp_data<-tmp_data[tmp_data$T==year,]
 
 
 n = 1000
 jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),bias=1.5) 
  lon_west=-17
  lon_east=54
  lat_south=-36.6
  lat_north=38
  par(bg="white")
  
  
  plot(1,1,type='n',xlim=c(lon_west,lon_east),ylim=c(lat_south,lat_north),xlab='',ylab='',axes=F,main=title)

  u <- par("usr")
  
  rect(u[1], u[3], u[2], u[4], col = rgb(180,225,255,maxColorValue=255), border = "black") 
  plot(World,add=T,col='white')

#### plot points
  my_cols<-match.cols(tmp_data$P,1000,jet.colors)
  
  points(tmp_data$X,tmp_data$Y,col=my_cols,pch=16,cex=0.7)
  
  
  plot(Africa,add=T)
  plot(Water,col = rgb(180,225,255,maxColorValue=255),add=T) 
  r.range <- c(0,1)
  leg.locs<-signif(seq(r.range[1],r.range[2],length.out=100),2)
  cols<-match.cols(leg.locs,n,jet.colors)
  label<-c(as.character(signif(r.range[1],2)),as.character(signif(r.range[2],2)))
  library(SDMTools)
  # xmin ymin xmax ymax
  rect(-17, -13, -8.5, -37, col = "white", border = "black") 
  pnts = cbind(x =c(-16.5,-12.5,-12.5,-16.5), y =c(-15,-15,-35.6,-35.6))
  legend.gradient(pnts,cols,label,title="")
  box()
}


relative_gain<-function(use,access){
	#use
	usegap=(access-use)/access # (relative proportion population with access to an ITN who do not sleep under it, calculated as (ACCESS-USE)/ACCESS )
	usegap[usegap>0]=0 #closing the use gap
	closeduserate=1-usegap # (proportion population with access to an ITN who do sleep under it, calculated as 1-USEGAP)
	### calculation of use rate as 1-(access-use)/access = use/access
	userate=use/access

	#access gain shows that under the current use rate, if access was 1 everywhere was is the level of use observed
	accessgain = userate #percentage point increase in coverage that could be achieved with universal access, but under existing USERATE):
	accessgain[accessgain>=1]=1

	# use rate gain shows that under current levels of access, if we closed the use gap what would the level of use be
	userategain = closeduserate*access #(percentage point increase in coverage that could be achieved with 100% use of existing nets, but under existing ACCESS):

	#relativegain=userategain/accessgain # (relative potential gain of increasing use rate alone versus increasing access alone) 
	relativegain=(accessgain-userategain)/accessgain # (relative potential gain of increasing use rate alone versus increasing access alone) 

	return(relativegain)
}

rel_gain_hist<-function(){
	h<-c()
	P=raster::stack(paste0('/home/backup/ITNcube/ITN_',2000:2016,'.ACC.tif'))
	P5=stack(paste0('/home/backup/ITNcube/ITN_',2000:2016,'.USE.tif'))
	ranges<-matrix(nrow=nlayers(P5),ncol=2)
	for(i in 1:nlayers(P5)){
		tmp<-relative_gain(P5[[i]],P[[i]])
		tmp<-getValues(tmp)
		tmp<-tmp[!is.na(tmp)]
		h<-c(h,tmp)
	}
	return(h)
	
}

rel_gain_range<-function(){

	P=raster::stack(paste0('/home/backup/ITNcube/ITN_',2000:2016,'.ACC.tif'))

	P5=stack(paste0('/home/backup/ITNcube/ITN_',2000:2016,'.USE.tif'))
	ranges<-matrix(nrow=nlayers(P5),ncol=2)
	for(i in 1:nlayers(P5)){
		tmp<-relative_gain(P5[[i]],P[[i]])
		tmp<-range(getValues(tmp),na.rm=T)
		ranges[i,]<-tmp
	}
	return(c(min(ranges[,1]),max(ranges[,2])))
	
}

dev_range<-function(){
	P=raster::stack(paste0('/home/backup/ITNcube/ITN_',2000:2016,'.DEV.tif'))
	ranges<-matrix(nrow=nlayers(P),ncol=2)
	for(i in 1:nlayers(P)){
		ranges[i,]<-ranges[i,]<-cbind(cellStats(P[[i]],stat='min'),cellStats(P[[i]],stat='max'))
	}
	return(c(min(ranges[,1]),max(ranges[,2])))

}

gap_range<-function(){
	P=raster::stack(paste0('/home/backup/ITNcube/ITN_',2000:2016,'.GAP.tif'))
	ranges<-matrix(nrow=nlayers(P),ncol=2)
	for(i in 1:nlayers(P)){
		ranges[i,]<-ranges[i,]<-cbind(cellStats(P[[i]],stat='min'),cellStats(P[[i]],stat='max'))
	}
	return(c(min(ranges[,1]),max(ranges[,2])))

}
#gainh<-rel_gain_hist() min value -2e1 max=1
gainr<-c(-2,1)
#gainr<-rel_gain_range()
#devr<-dev_range()
#gapr<-gap_range()

data<-read.csv('/home/backup/ITNcube/ITN_final_clean_access_11thSept2017.csv')



file.remove('/home/backup/ITNcube/access.pdf')
filename<-'/home/backup/ITNcube/access.pdf'
pdf(filename,width=18.61, height=10.34,onefile=T) 

for(i in time_points){
	P=raster(paste0('/home/backup/ITNcube/ITN_',i,'.ACC.tif'),NAflag=-9999)
#	P2=raster(paste0('/home/backup/ITNcube/ITN_',i,'.DEV.tif'),NAflag=-9999)
	P3=raster(paste0('/home/backup/ITNcube/ITN_',i,'.MEAN.tif'),NAflag=-9999)
#	P4=raster(paste0('/home/backup/ITNcube/ITN_',i,'.GAP.tif'),NAflag=-9999)
	P5=raster(paste0('/home/backup/ITNcube/ITN_',i,'.USE.tif'),NAflag=-9999)
	relg<-relative_gain(P5,P)
 	relg_decision<-relg
    relg_decision[relg_decision>0]=1
	relg_decision[relg_decision<0]=-1
	relg_decision[relg_decision==0]=0
	
	P[1:5]=1;P[6:10]=0
	P3[1:5]=1;P3[6:10]=0
	P5[1:5]=1;P5[6:10]=0
#	P2[1:5]=devr[1];P2[6:10]=devr[2]
#	P4[1:5]=gapr[1];P4[6:10]=gapr[2]
	relg[1:5]=gainr[1];relg[6:10]=gainr[2]
	relg[relg<gainr[1]]=gainr[1]
	P6<-(P-P5)/P
	P6=1-P6
	P6[P6<1e-2]=0
	
	par(mfrow=c(2,4),mar=c(1,1,1,1))
	Plot_raster(P,1,"ITN access")
#	Plot_raster(P2,1,"ITN access deviation")
	Plot_raster(P3,1,"ITN access mean")
#	Plot_raster(P4,1,"ITN use gap")
	Plot_raster(P5,1,"ITN use")
	Plot_raster(P6,1,"Off Acces")
	Plot_raster(relg,0.3,"ITN relative gain")
	Plot_raster(relg_decision,1,"ITN decision gain")
	Plot_points(data,i,'use','Raw use data')
	Plot_points(data,i,'access','Raw access data')
	

}
dev.off()

#use
use=P3
access=P	
r=use-access
r[r<=0]=0
l=r+P


#access
use=P3
access=P	
r=use-access
l2=1+r
l2[l2>=1]=1



#use
use=P3 #(proportion population sleeping under a net)
access=P #(proportion population with access to an ITN in their household)	 
usegap=(access-use)/access # (relative proportion population with access to an ITN who do not sleep under it, calculated as (ACCESS-USE)/ACCESS )
usegap[usegap>0]=0 #closing the use gap
closeduserate=1-usegap # (proportion population with access to an ITN who do sleep under it, calculated as 1-USEGAP)
### calculation of use rate as 1-(access-use)/access = use/access
userate=use/access

#access gain shows that under the current use rate, if access was 1 everywhere was is the level of use observed
accessgain = userate #percentage point increase in coverage that could be achieved with universal access, but under existing USERATE):
accessgain[accessgain>=1]=1

# use rate gain shows that under current levels of access, if we closed the use gap what would the level of use be
userategain = closeduserate*access #(percentage point increase in coverage that could be achieved with 100% use of existing nets, but under existing ACCESS):

#relativegain=userategain/accessgain # (relative potential gain of increasing use rate alone versus increasing access alone) 
relativegain=(accessgain-userategain)/accessgain # (relative potential gain of increasing use rate alone versus increasing access alone) 





par(mfrow=c(1,2))

plot(relativegain,col=jet.colors(1000),main='relative gain')

 relativegain[relativegain>0]=1
 relativegain[relativegain<0]=-1
 relativegain[relativegain==0]=0

plot(relativegain,main='relative gain decision')
 



get.cn.estimates<-function(r,pop){
	database<-read.csv('/home/backup/country_table.csv')
	cn<-raster('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif') # master country layer
	NAvalue(cn)=-9999
	try(if(extent(cn)!=extent(r)) stop("extents do not match"))
	limits<-raster('/home/drive/cubes/5km/Pf_limits/Pf_limits.tif')
	NAvalue(limits)=-9999
	# stable only limits
	limits[limits==2]=2
	limits[limits==9999 | limits==0 | limits==1]=NA
	# no na values in county population or limits
	NAsum<-cn+pop+limits+r
	NAsum<-getValues(NAsum)
	w=!is.na(NAsum)	
	
	mat<-rep(NA,nrow=nrow(database))
	# ger values
	cn<-getValues(cn)
	cn=cn[w]
		
	pop<-getValues(pop)
	pop=pop[w]
	
	r=getValues(r)
	r=r[w]

		
	for(i in 1:nrow(database)){
		gaul<-database[i,'GAUL_CODE']
		cpop<-pop[cn==gaul]
		popsum<-sum(cpop,na.rm=TRUE)
		cpop<-cpop/popsum
		
		mat[i]=weighted.mean(r[cn==gaul],w=cpop,na.rm=T)


	}
	return(mat)	
}

	
library(doParallel)
registerDoParallel(62)
tables<-foreach(i=time_points) %dopar% {	
	P=raster(paste0('/home/backup/ITNcube/ITN_',i,'.ACC.tif'),NAflag=-9999)
	P2=raster(paste0('/home/backup/ITNcube/ITN_',i,'.MEAN.tif'),NAflag=-9999)
	P3=raster(paste0('/home/backup/ITNcube/ITN_',i,'.USE.tif'),NAflag=-9999)
	pop<-raster(paste0('/home/drive/cubes/5km/AfriPop/',i,'.total.population.tif'),NAflag=-9999)
	P4=(P-P2)/P2 # relative access deviation
	P5=(P3-P)/P # relative use gap
	P6=(P-P2) # absolute access deviation
	P7=(P3-P) # absolute use gap
	
	t1<-get.cn.estimates(P,pop)
	t2<-get.cn.estimates(P2,pop)
	t3<-get.cn.estimates(P3,pop)
	t4<-get.cn.estimates(P4,pop)
	t5<-get.cn.estimates(P5,pop)
	t6<-get.cn.estimates(P6,pop)
	t7<-get.cn.estimates(P7,pop)
	
	return(cbind(t1,t2,t3,t4,t5,t6,t7))
	
}	

#####use as a function of access

file.remove('/home/backup/ITNcube/Summary.pdf')
filename<-'/home/backup/ITNcube/Summary.pdf'
pdf(filename,width=11.59, height=10.06,onefile=T) 
	i=2013
for(i in time_points){

	P=raster(paste0('/home/backup/ITNcube/ITN_',i,'.ACC.tif'),NAflag=-9999)
	P2=raster(paste0('/home/backup/ITNcube/ITN_',i,'.MEAN.tif'),NAflag=-9999)
	P5=raster(paste0('/home/backup/ITNcube/ITN_',i,'.USE.tif'),NAflag=-9999)
	
	
	P[1:5]=1;P[6:10]=0
	P2[1:5]=1;P2[6:10]=0
	P5[1:5]=1;P5[6:10]=0

	P6<-(P-P5)/P
	P6=1-P6
	P6[P6<1e-2]=0

	par(mfrow=c(2,2),mar=c(1,1,1,1))
	Plot_raster(P,1,"ITN access")
	Plot_raster(P2,1,"ITN Stock and flow mean")
	Plot_raster(P5,1,"ITN use")
	Plot_raster(P6,1,"Use as function of Access")
	
}
	dev.off()



# this routine adjust use numbers to zero in bad years and label files the same as before for backward compatibility
indicators<-read.csv('/home/backup/ITNcube/indicators_access_qtr_new.csv')
indicatorsr<-read.csv('/home/backup/ITNcube/indicators_access_qtr_new.csv')

ind<-matrix(nrow=nrow(indicators),ncol=17)
for(i in 1:nrow(indicators)){
	ind[i,]<-tapply(as.numeric(indicators[i,2:ncol(indicators)]),c(rep(1:(72/4), each = 4),18),mean)[1:17]
}
rownames(ind)<-indicators[,1]
colnames(ind)<-seq(2000,2016,by=1)
indicators<-ind

POPULATIONS<-read.csv('/home/backup/ITNcube/country_table_populations.csv') # load table to match gaul codes to country names
names<-as.character(rownames(ind))
#for(i in 1:nrow(indicators)){
#	names[i]=as.character(POPULATIONS[as.character(POPULATIONS$NAME)==names[i],'COUNTRY_ID']) # get 3 letter country codes
#}
rownames(ind)<-as.character(names)
threshold=0.02
indmat<-matrix(data=0,nrow=nrow(indicators),ncol=ncol(indicators))
for(i in 1:nrow(indicators)){
	for(j in 1:ncol(indicators)){
	 	if(indicators[i,j]<threshold){
	  		indmat[i,j]=1
	  	}
		if(indicators[i,j]>threshold){
			break
		}
	}
}
colnames(indmat)<-seq(2000,2016,by=1)
rownames(indmat)<-names

times<-seq(2000,2016,by=1)
st<-stack(paste('/home/backup/ITNcube/ITN_',times,'.USE.tif',sep=""))
stv<-getValues(st)
cn<-raster(paste('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif',sep=""))
cnv<-getValues(cn)
for(i in 1:nrow(indmat)){
	gaul<-POPULATIONS[as.character(POPULATIONS$NAME)==rownames(ind)[i],'GAUL_CODE']
	tmp<-cn
	wh<-colnames(indmat)[indmat[i,]==1]
	wh<-as.numeric(wh)-1999
	stv[cnv==gaul,wh]<-0
}

for(i in 1:nlayers(st)){
	st[[i]]<-setValues(st[[i]],stv[,i])
}

for(i in 1:nlayers(st)){
	writeRaster(st[[i]],paste('/home/drive/cubes/5km/ITN/',times[i],'.ITN.use.yearavg.new.adj.tif',sep=""),NAflag=-9999,overwrite=TRUE)

}

# create zero years
st0<-st[[1]]
NAvalue(st0)<--9999
st0[!is.na(st0)]=0
writeRaster(st0,paste('/home/drive/cubes/5km/ITN/',1999,'.ITN.use.yearavg.new.adj.tif',sep=""),NAflag=-9999,overwrite=TRUE)
writeRaster(st0,paste('/home/drive/cubes/5km/ITN/',1998,'.ITN.use.yearavg.new.adj.tif',sep=""),NAflag=-9999,overwrite=TRUE)
writeRaster(st0,paste('/home/drive/cubes/5km/ITN/',1997,'.ITN.use.yearavg.new.adj.tif',sep=""),NAflag=-9999,overwrite=TRUE)
writeRaster(st0,paste('/home/drive/cubes/5km/ITN/',1996,'.ITN.use.yearavg.new.adj.tif',sep=""),NAflag=-9999,overwrite=TRUE)
writeRaster(st0,paste('/home/drive/cubes/5km/ITN/',1995,'.ITN.use.yearavg.new.adj.tif',sep=""),NAflag=-9999,overwrite=TRUE)

