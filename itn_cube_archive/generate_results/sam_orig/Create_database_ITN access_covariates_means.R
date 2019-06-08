library(zoo)
library(raster)
library(VGAM)
load('/home/backup/ITNcube/prop0prop1.rData')
emplogit<-function(Y,N) {top=Y*N+0.5;bottom=N*(1-Y)+0.5;return(log(top/bottom))} # empirical logit number
emplogit2<-function(Y,N) {top=Y+0.5;bottom=N-Y+0.5;return(log(top/bottom))} # empirical logit two numbers
######################################### FUNCTIONS ###########################################################
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
				st<-paste('<',points$Survey[x],'-',points$year[x],'>')
				print(paste('-->Point',st ,'has been repositioned at distance of',min(dist)/1000,"km","at(",points[x,'lon'],",",points[x,'lat'],")"))  	
				corindex<-c(corindex,x)
		  } else {
				st<-paste('<',points$Survey[x],'-',points$year[x],'>')
			print(paste("-->Point ",st,points[x,'lon'],points[x,'lat'],"is not positioned properly"))
		  }
		}
	}
	return(points)

	}
	

aggregate.data<-function(f,cn){
	# for spatial only data snap points in 5km pixel to centroid and then aggregate number examined and positive
	r<-cn
	ce<-cellFromXY(r,f[,c('lon','lat')])  # get cell numbers
	f$cellnumber<-ce # update the cell numbers
	cexy<-xyFromCell(r,ce)	# get xy from cells
	catData<-rep(NA,length(ce))		
	# concatenate cellX-cellY-month-year for unique entry searching
	library(doParallel)
	registerDoParallel(60)	
	catData<-foreach(i=1:length(ce)) %dopar% {	
		 paste(cexy[i,1],cexy[i,2],f$yearqtr[i],sep=':')
	}
	catData<-unlist(catData)
	un<-unique(catData)		#find points with unique lat,lon and time
	#### required columns of f
	fcolapse<-matrix(nrow=length(un),ncol=ncol(f))	#create colapse matrix to sum aggregated data
	colnames(fcolapse)<-colnames(f)				#label		
	fcolapse<-as.data.frame(fcolapse)	
	f[,'lon']<-cexy[,1]	# add fixed lons
	f[,'lat']<-cexy[,2]	#add fixed lats
	library(doParallel)
	registerDoParallel(60)	
	aggregate<-foreach(i=1:length(un), .combine=rbind) %dopar% {	
		tmp<-f[catData==un[i],] #get entries (repeated or single)
		ret<-tmp[1,] # get first entry. appart from + and tot and date all entries should be the same
		if(nrow(tmp)==1){
			# do nothing
		} else {
			ret['P']<-as.integer(sum(tmp[,'P'])) #add the positives
			ret['N']<-as.integer(sum(tmp[,'N'])) # add the examined	 
		}
		return(ret)	
	}
	return(aggregate)
}


calc_access_matrix<-function(accdate1,func0,func1,hh,nc=10,nr=40){
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

		allind<-get.access.mat(matModel,nc,nr)	
		
		acc<-rep(NA,10)
		for(i in 1:10){
			acc[i]<-sum(allind[[1]][,i])/sum(allind[[2]][,i])
		}
		acc<-c(acc,sum(allind[[1]])/sum(allind[[2]]))
		return(acc)
}

get.access.mat<-function(mat,nc,nr){
	cnames<-colnames(mat)<-1:nc
	rnames<-rownames(mat)<-0:nr	

	ind3mat<-matrix(data=0,nrow=nrow(mat),ncol=ncol(mat))
	ind3tot<-matrix(data=0,nrow=nrow(mat),ncol=ncol(mat))

	for(x in 1:ncol(mat)){
		for(y in 1:nrow(mat)){
			ind3mat[y,x]<-2*(y-1)*mat[y,x] 
			ind3tot[y,x]<-x*mat[y,x] # total people
			if(ind3mat[y,x]>ind3tot[y,x]) ind3mat[y,x]=ind3tot[y,x] #cap on estimate 
		}
	}
	l<-list()
	l[[1]]=ind3mat
	l[[2]]=ind3tot
	return(l)
}

KEY=read.csv('/home/backup/ITNcube/KEY_080817.csv')
#
#print(paste('***OUTPUT: loading DHS surveys - both household and cluster'))
#HH<-read.csv('/home/backup/ITNcube/DHS_MIS_all_080817.csv',stringsAsFactors=F)
#print(paste('***OUTPUT:Link household and cluster to calculate access'))
#Surveys<-as.character(unique(HH$Survey.hh))
#
##### add some of the extra data bonnie got like zambia 2010
#HH1<-read.csv('/home/backup/ITNcube/Net details aggregated by cluster combined19March2015_all.csv',stringsAsFactors=F)
#HH2<-read.csv('/home/backup/ITNcube/Net details aggregated by household combined19March2015_all.csv',stringsAsFactors=F)
#
#harrylist<-as.character(KEY[KEY[,1]%in%Surveys,2])
#bonnielist<-as.character(unique(HH2$Survey.hh))
#additionalSVY<-bonnielist[!bonnielist%in%harrylist] # these survey are in bonnies files but not harry
#additionalSVY<-additionalSVY[!additionalSVY%in%"Zambia 2010"]
#additionalSVY<-additionalSVY[!additionalSVY%in%"Zambia 2012"]
#
#
#
#HH1<-HH1[HH1$Survey%in%additionalSVY,]
#HH2<-HH2[HH2$Survey.hh%in%additionalSVY,] #restrict to aditional surveys
#HH2$hhid=NA
#HH2$latitude=NA
#HH2$longitude=NA
#HH2<-HH2[,colnames(HH)]
#
#for(i in 1:length(additionalSVY)){
#	tmp<-HH1[HH1$Survey==additionalSVY[i],]
#	un<-unique(tmp$Cluster.number)
#	tmp2<-HH2[HH2$Survey.hh==additionalSVY[i],]
#
#	for(j in 1:length(un)){
#		tmp2$longitude[tmp2$Cluster.hh==un[j]]=as.numeric(tmp$Long.cluster[tmp$Cluster.number==un[j]])
#		tmp2$latitude[tmp2$Cluster.hh==un[j]]=as.numeric(tmp$Lat.cluster[tmp$Cluster.number==un[j]])
#	}
#	HH<-rbind(HH,tmp2) # add bonnies shitty data
#}
#
#
#HH$year<-as.numeric(as.yearmon(paste(HH[, "year"], "-", HH[, "month"],sep=""))) # calculate numeric date
#################################### for DHS
#HH<-HH[!is.na(HH$n.individuals.that.slept.in.surveyed.hhs),]
#HH<-HH[!is.na(HH$n.ITN.per.hh),]
#HH<-HH[!is.na(HH$n.individuals.that.slept.under.ITN),]
#HH<-HH[!is.na(HH$latitude),]
#HH<-HH[!is.na(HH$longitude),]
#rm(HH1)
#rm(HH2)
Country.list<-Cout
HH<-read.csv('/home/backup/ITNcube/ALL_HH_Data_20112017.csv',stringsAsFactors=F)

HH<-HH[HH$ISO3%in%Country.list,]



zerohh<-is.na(HH$n.individuals.that.slept.in.surveyed.hhs)
HH<-HH[!zerohh,] # remove entries with NA household size

zerohh<-(HH$n.individuals.that.slept.in.surveyed.hhs)==0
HH<-HH[!zerohh,] # remove entries with zero household size

Surveys<-as.character(unique(HH$Survey.hh)) #update survey vector
Surveys<-Surveys[!is.na(Surveys)]

times<-seq(0.0,17,0.25)+2000
a=HH$n.individuals.that.slept.in.surveyed.hhs # number in household
b=HH$n.ITN.per.hh # number of nets
d=b*2 # number covered by nets
d[d>a]=a[d>a] # bounding to not exceed the number in the house
e=HH$n.individuals.that.slept.under.ITN # use count


library(doParallel)
registerDoParallel(62)		
output<-foreach(i=1:length(Surveys),.combine=rbind) %dopar% { #survey loop
	svy<-Surveys[i] # store survey
	
	cc<-HH[HH$Survey.hh==svy,'ISO3'] # get country
	p<-out[[which(Country.list==cc[1])]]
	p0<-p[,,1] # get p0 parameters
	p1<-p[,,2] # get p1 parameters	
	
	func0<-list()
	func1<-list()
	# make these into continuous functions
	for(f in 1:ncol(p0)){
		eval(parse(text=paste0('func0[[',f,']]<-approxfun(times,p0[,',f,'])')))
	}
	for(f in 1:ncol(p0)){
		eval(parse(text=paste0('func1[[',f,']]<-approxfun(times,p1[,',f,'])')))
	}

	tmp_hh=HH[HH$Survey.hh==svy,] # subset hosehold 
	un<-unique(tmp_hh$Cluster.hh) # get unique cluster IDs for that household survey
	### house hold size as measured by the number of people sleeping under a net the night before
	hh_val<-rep(NA,100) 
	tmp_hh$sample.w<-tmp_hh$sample.w/sum(tmp_hh$sample.w) # get sample weights
	for(xx in 1:100){
		hh_val[xx]<-sum(tmp_hh[tmp_hh$n.individuals.that.slept.in.surveyed.hhs==xx,'sample.w'])	# proportion in a given household size	
	}
	hh_size<-hh_val[1:10] # 1:10
	hh_size[10]<-sum(hh_val[10:length(hh_val)])	# >10
	hh<-as.vector(as.numeric(hh_size))
	names(hh)<-1:10	
	#store matrix
	fields<-c("Survey" ,"lat" , "lon",'year','P','N','Pu','T','Amean','Tmean','gap','gap2','gap3')
	tmp<-matrix(NA,nrow=length(un),ncol=length(fields))
	colnames(tmp)<-fields
	tmp<-as.data.frame(tmp)

	for(j in 1:length(un)){ #for each cluster
		wh<-HH$Survey.hh==svy & HH$Cluster.hh==un[j] # get the right HH values



		tmp[j,'Survey']<-svy
		tmp[j,'lat']<-mean(HH[wh,'latitude'],na.rm=TRUE)
		tmp[j,'lon']<-mean(HH[wh,'longitude'],na.rm=TRUE)

		tmp[j,'P']<-sum(d[wh]) #access cluster count
		tmp[j,'N']<-sum(a[wh]) #number in cluster count
		tmp[j,'year']<-mean(HH$year[wh])
		tmp[j,'Pu']<-sum(e[wh]) #number in cluster count
		tmp[j,'T']<-sum(b[wh]) #number in cluster count
		tmp[j,'gap']<-mean(((tmp[j,'P']/tmp[j,'N'])-(tmp[j,'Pu']/tmp[j,'N']))/(tmp[j,'P']/tmp[j,'N'])) # (access-use)/access
		tmp[j,'gap2']<-mean(emplogit2(tmp[j,'P'],tmp[j,'N'])-emplogit2(tmp[j,'Pu'],tmp[j,'N'])) # emplogit difference of access-use
		
		tmp2<-HH[wh,]
		
		acc<-acc_mean<-agap<-ugap<-rep(NA,sum(wh))
		P<-d[wh]
		Pu<-e[wh]
		N<-a[wh]
		for(k in 1:nrow(tmp2)){
			y<-tmp2$year;h<-tmp2$n.individuals.that.slept.in.surveyed.hhs;
			h[h>10]=10	
			access<-calc_access_matrix(y[k],func0,func1,hh)
			acc[k]<-access[h[k]]
			acc_mean[k]<-access[11]		
			ugap[k]<- 1-(Pu[k]/P[k])
		}
		tmp[j,'Amean']<-mean(acc,na.rm=T)
		tmp[j,'Tmean']<-mean(acc_mean,na.rm=T)	
		tmp[j,'gap3']<-mean(ugap,na.rm=T)

	}
	return(tmp)
}

HH1<-output

data<-HH1
#print(paste('**OUTPUT MESSAGE** remove points with no cooridnates: there are - ',nrow(data[!complete.cases(data),])))
data<-data[complete.cases(data),]

print(paste('**OUTPUT MESSAGE** remove points with 0 lat 0 lon: there are - ',nrow(data[data$lat==0 & data$lon==0,])))
data<-data[data$lat!=0 & data$lon!=0,]

# Check for invalid points, and attempt to reposition them
cn<-raster('/home/drive/cubes/5km/Admin/african_cn5km_2013_no_disputes.tif') # master country layer
NAvalue(cn)=-9999
data$yearqtr<-as.numeric(as.yearqtr(data$year))


print(paste('**OUTPUT MESSAGE** Atempting to reposition points'))
data<-reposition.points(cn,data,4)


print(paste('-->Total number of household points', nrow(HH)))
print(paste('-->Total number of cluster points', nrow(data)))
print(paste('-->Total number countries', length(unique(HH$Country))))
print(paste('-->Total number of surveys', length(unique(HH$Survey.hh))))



print(paste('**OUTPUT MESSAGE** Aggregating data at same pixel-quarter'))
data<-aggregate.data(data,cn)

print(paste('**OUTPUT MESSAGE** get floored year '))

data$flooryear<-floor(data$year)



write.csv(data,'/home/backup/ITNcube/ITN_final_clean_access_20thNov2017.csv',row.names=FALSE)
