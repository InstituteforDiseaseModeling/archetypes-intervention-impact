
## straight from sam (/home/backup/Space Time Malaria/) used (?) in acc_deviation_refactored.r
## todo: add to acc_dev_functions

########################################################################################################################  
##### functions #######  ########################################################################################## 
##create.list.cov<-function(data){
#	est<-c( list(upr_u=data[,"upr_u"]),
#	list(upr_p=data[,"upr_p"]),
#	list(prec57a0=data[,"prec57a0"]),
#	list(prec57a1=data[,"prec57a1"]),
#	list(prec57a2=data[,"prec57a2"]),
#	list(prec57mx=data[,"prec57mx"]),
#	list(prec57p1=data[,"prec57p1"]),
#	list(prec57p2=data[,"prec57p2"]),
#	list(tempaucpf=data[,"tempaucpf"]),
#	list(wd0103a0=data[,"wd0103a0"]),
#	list(wd0103mx=data[,"wd0103mx"]),
#	list(wd0107a1=data[,"wd0107a1"]),
#	list(wd0107a3=data[,"wd0107a3"]),
#	list(wd0107mn=data[,"wd0107mn"]),
#	list(wd0107p1=data[,"wd0107p1"]),
#	list(wd0114a0=data[,"wd0114a0"]),
#	list(wd0114a3=data[,"wd0114a3"]),
#	list(wd0114mx=data[,"wd0114mx"]),
#	list(wd0114p1=data[,"wd0114p1"]),
#	list(wd0114p2=data[,"wd0114p2"]),
#	list(time=data[,"months"]))
#	return(est)
#}
#
summary.obsmat<-function(A){
	dimensions<-dim(A)
	loc<-table(rowSums(A>0))
	unused<-table(colSums(A)>0)

	print(paste("A matrix has", dimensions[1], "rows (points) and",dimensions[2],"columns (nodes)"))
	print(paste(loc[1], "points lie on a vertex and", loc[2], "points lie in a triangle"))
	print(paste("There are",unused[1], "unused nodes and", unused[2],"used nodes"))

}


create.list.cov<-function(data){
	names<-colnames(data)
	remove<-c('x','y','z','response','examined','prevalence')
	names<-names[! names %in% remove]

	data<-data[,names]
	est<-as.list(data)
	return(est)
}


extract.evi<-function(d,lag){
 	numericmon<-c('04','08','12','02','01','07','06','03','05','11','10','09')
    numericyear<-as.character(2000:2012)
	evi<-rep(NA,nrow(d))
	response<-rep(NA,nrow(d))
	N<-rep(NA,nrow(d))
	Npos<-rep(NA,nrow(d))
	if(lag==0){
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						lm=m;ly=y
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$examined)
						Npos[tmp$id]<-as.integer(tmp$positive)
						print(paste(m,y,'with no lag'))
					}
				}
				j=j+1
			}
			i=i+1
		}		
	} else if(lag==1){
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)
						}else {
						  lm=sortedmon[j-1]			
						  ly=y
						}
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$examined)
						Npos[tmp$id]<-as.integer(tmp$positive)				
						print(paste(m,y,'has been lagged by 1 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}	
	} else if(lag==2){
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='11'		
						   ly=as.character(as.numeric(y)-1)
						} else if(m=='02'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)				
						}else {
						  lm=sortedmon[j-2]			
						  ly=y
						}
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$examined)
						Npos[tmp$id]<-as.integer(tmp$positive)					
						print(paste(m,y,'has been lagged by 2 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}	
	} else if(lag==3){
		# EVI new lagged by 3
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='10'		
						   ly=as.character(as.numeric(y)-1)
						} else if(m=='02'){
						   lm='11'		
						   ly=as.character(as.numeric(y)-1)				
						} else if(m=='03'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)				
						}else {
						  lm=sortedmon[j-3]			
						  ly=y
						}
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$examined)
						Npos[tmp$id]<-as.integer(tmp$positive)						
						print(paste(m,y,'has been lagged by 3 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}
	} else if(lag==4){
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='09'		
						   ly=as.character(as.numeric(y)-1)
						} else if(m=='02'){
						   lm='10'		
						   ly=as.character(as.numeric(y)-1)				
						} else if(m=='03'){
						   lm='11'		
						   ly=as.character(as.numeric(y)-1)				
						}else if(m=='04'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)				
						}else {
						  lm=sortedmon[j-4]			
						  ly=y
						}
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$examined)
						Npos[tmp$id]<-as.integer(tmp$positive)					
						print(paste(m,y,'has been lagged by 4 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}	
	
	} else if(lag==5){
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='08'		
						   ly=as.character(as.numeric(y)-1)
						} else if(m=='02'){
						   lm='09'		
						   ly=as.character(as.numeric(y)-1)				
						} else if(m=='03'){
						   lm='10'		
						   ly=as.character(as.numeric(y)-1)				
						}else if(m=='04'){
						   lm='11'		
						   ly=as.character(as.numeric(y)-1)				
						}else if(m=='05'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)				
						}else {
						  lm=sortedmon[j-5]			
						  ly=y
						}
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$positive)
						Npos[tmp$id]<-as.integer(tmp$examined)						
						print(paste(m,y,'has been lagged by 5 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}		
	} else if(lag==6){
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='07'		
						   ly=as.character(as.numeric(y)-1)
						} else if(m=='02'){
						   lm='08'		
						   ly=as.character(as.numeric(y)-1)				
						} else if(m=='03'){
						   lm='09'		
						   ly=as.character(as.numeric(y)-1)				
						}else if(m=='04'){
						   lm='10'		
						   ly=as.character(as.numeric(y)-1)				
						}else if(m=='05'){
						   lm='11'		
						   ly=as.character(as.numeric(y)-1)				
						}else if(m=='06'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)				
						}else {
						  lm=sortedmon[j-6]			
						  ly=y
						}
						r<-raster(paste('/home/drive/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						N[tmp$id]<-as.integer(tmp$positive)
						Npos[tmp$id]<-as.integer(tmp$examined)						
						print(paste(m,y,'has been lagged by 6 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}	
	}
	return(cbind(evi,response,N,Npos))
}

get.pred.locs<-function(wh){

	if(is.numeric(wh)){
		r<-raster('/home/likewise-open/ZOO/bras2280/Space Time Malaria/mcmcpfreg5k_CLEAN.flt')
		r[r!=wh]<-NA
	} else {
		r<-raster('/home/likewise-open/ZOO/bras2280/Space Time Malaria/PF_LIMITS.tif')
	}
	pred_val<-getValues(r)
	w<-is.na(pred_val)
	index<-1:length(w)
	index<-index[!w]
	pred_locs<-xyFromCell(r,1:ncell(r))  
	pred_locs<-pred_locs[!w,] 

	colnames(pred_locs)<-c('longitude','latitude')
	pred_locs<-ll.to.xyz(pred_locs)
	return(cbind(pred_locs))
}

get.pred.data<-function(wh){
	raster.fnames<-c("/home/drive/mastergrids/cleandata/tempsuitability/pf/tempaucpf_5k_CLEAN.flt ",
	"/home/drive/mastergrids/cleandata/urbanareas/grumpALPHA/upr_u_5km.flt",
	"/home/drive/mastergrids/cleandata/urbanareas/grumpALPHA/upr_p_5km.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0114a0_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/tfaWCLIMprec/prec57a0_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/tfaWCLIMprec/prec57a1_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/tfaWCLIMprec/prec57a2_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/tfaWCLIMprec/prec57mx_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/tfaWCLIMprec/prec57p1_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/tfaWCLIMprec/prec57p2_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0103a0_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0103mx_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0107a1_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0107a3_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0107p1_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0107mn_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0114mx_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0114p1_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0114p2_5k_CLEAN.flt",
	"/home/drive/mastergrids/cleandata/AIPavhrr1km/wd0114a3_5k_CLEAN.flt")
	fnames<-c("tempaucpf",
					"upr_u",
					"upr_p",
					"wd0114a0",
					"prec57a0",
					"prec57a1",
					"prec57a2",
					"prec57mx",
					"prec57p1",
					"prec57p2",
					"wd0103a0",
					"wd0103mx",
					"wd0107a1",
					"wd0107a3",
					"wd0107p1",
					"wd0107mn",
					"wd0114mx",
					"wd0114p1",
					"wd0114p2",
					"wd0114a3")				
	s<-stack(raster.fnames)
	if(is.numeric(wh)){
		r<-raster('/home/likewise-open/ZOO/bras2280/Space Time Malaria/mcmcpfreg5k_CLEAN.flt')
		r[r!=wh]<-NA
	} else {
		r<-raster('/home/likewise-open/ZOO/bras2280/Space Time Malaria/PF_LIMITS.tif')
	}
	pred_val<-getValues(r)
	w<-is.na(pred_val)
	index<-1:length(w)
	index<-index[!w]
	pred_locs<-xyFromCell(r,1:ncell(r))  
	pred_locs<-pred_locs[!w,] 
	pred_data<-matrix(nrow=length(index),ncol=length(raster.fnames))
	for(i in 1:length(raster.fnames)){
		print(paste("On raster stack number",i))
		pred_data[,i]<-getValues(s[[i]])[index]
	}
#	pred_data<-extract(s,index)
	colnames(pred_data)<-fnames
	colnames(pred_locs)<-c('longitude','latitude')
	pred_locs<-ll.to.xyz(pred_locs)
	return(cbind(pred_locs,pred_data))
}

get.pred.data.generic<-function(wh,raster.fnames,fnames){
	
	s<-stack(raster.fnames)
	if(is.numeric(wh)){
		r<-raster('/home/likewise-open/ZOO/bras2280/Space Time Malaria/mcmcpfreg5k_CLEAN.flt')
		print(paste("Using MCMC region",wh))
		r[r!=wh]<-NA
		print(paste("Trimming raster to MCMC region"))
		r<-trim(r)
	} else {
		print(paste("using full PF limits"))
		r<-raster('/home/likewise-open/ZOO/bras2280/Space Time Malaria/PF_LIMITS.tif')
	}
	pred_val<-getValues(r)
	w<-is.na(pred_val)
	index<-1:length(w)
	index<-index[!w]
	pred_locs<-xyFromCell(r,1:ncell(r))  
	pred_locs<-pred_locs[!w,] 
	pred_data<-matrix(nrow=length(index),ncol=length(raster.fnames))
	for(i in 1:length(raster.fnames)){
		print(paste("On raster stack number",i,"which is:",fnames[i]))
		pred_data[,i]<-getValues(s[[i]])[index]
	}
#	pred_data<-extract(s,index)
	colnames(pred_data)<-fnames
	colnames(pred_locs)<-c('longitude','latitude')
	pred_locs<-ll.to.xyz(pred_locs)
	return(cbind(pred_locs,pred_data))
}

validation.stats<-function(mod.pred,data.val,data.est,actual_est,actual_val){
	#validation ############################################################################
	index.val = inla.stack.index(stack,"val")$data
	lp.mean.val = mod.pred$summary.linear.predictor[index.val,"mean"]
	lp.sd.val = mod.pred$summary.linear.predictor[index.val,"sd"]
	lp.025.val = mod.pred$summary.linear.predictor[index.val,"0.025quant"]
	lp.975.val = mod.pred$summary.linear.predictor[index.val,"0.975quant"]
	
	lp.mean.val<-inv.logit(lp.mean.val)
	lp.sd.val<-inv.logit(lp.sd.val)
	lp.025.val<-inv.logit(lp.025.val)
	lp.975.val<-inv.logit(lp.975.val)
	
	#estimation ###########################################################################
	index.est = inla.stack.index(stack,"est")$data
	lp.mean.est = mod.pred$summary.linear.predictor[index.est,"mean"]
	lp.sd.est = mod.pred$summary.linear.predictor[index.est,"sd"]
	lp.025.est = mod.pred$summary.linear.predictor[index.est,"0.025quant"]
	lp.975.est = mod.pred$summary.linear.predictor[index.est,"0.975quant"]
	
	lp.mean.est<-inv.logit(lp.mean.est)
	lp.sd.est<-inv.logit(lp.sd.est)
	lp.025.est<-inv.logit(lp.025.est)
	lp.975.est<-inv.logit(lp.975.est)



	stats<-list()
	mean_actual=mean(lp.mean.val)
	stats$r2 <- 1-(sum((lp.mean.val-actual_val)^2)/sum((actual_val-mean_actual)^2))
	stats$dic<-mod.pred$dic$dic 
	stats$mse<- mean((lp.mean.val-actual_val)^2)
	stats$me<-mean((lp.mean.val-actual_val))
	stats$mae<-mean(abs(lp.mean.val-actual_val))
	stats$est.cor<-cor(lp.mean.est,actual_est)
	stats$val.cor<-cor(lp.mean.val,actual_val)
	stats$est.sd<-sd(lp.mean.est-actual_est)
	stats$val.sd<-sd(lp.mean.val-actual_val)
	stats$mean.val<-lp.mean.val
	stats$actual.val<-actual_val
	return(stats)
}


quantile.coverage<-function(mod.pred,data.est){
	quants<-seq(0.01,0.99,0.01)
	actual_est<-data.est[,'response']/data.est[,'examined']
	index.est = inla.stack.index(stack,"est")$data
	mm<-matrix(ncol=length(quants),nrow=length(index.est))
	stats<-list()
	for(i in index.est){
		mm[i,]<-inla.qmarginal(quants,mod.pred$marginals.linear.predictor[[i]])
	}
	mm<-inv.logit(mm)
	
	count<-rep(0,length(quants))
	for(i in 1:length(quants)){
		count[i]<-sum(actual_est<mm[,i])
	}
	count<-count/nrow(mm)
	stats$q<-quants
	stats$coverage<-count
	return(stats)
}

model.coverage.val<-function(mod.pred,data.val){
	library(boot)

	
	#valimation ###########################################################################
	index.val = inla.stack.index(stack,"val")$data
	actual_val<-data.val[,'response']/data.val[,'examined']
	stats<-list()
	
	qr<-function(q){return(cbind(.5*(1.-q),1-(.5*(1.-q))))}
	quants<-c(qr(0.95),qr(0.85),qr(0.75),qr(0.65),qr(0.55),qr(0.45),qr(0.35),qr(0.25),qr(0.15),qr(0.05))
	mm<-matrix(ncol=length(quants),nrow=length(index.val))
	for(i in index.val){
		mm[i,]<-inla.qmarginal(quants,mod.pred$marginals.linear.predictor[[i]])
	}
	mm<-inv.logit(mm)
	
	vec<-c()
	i<-1
	while(i < (length(quants))){
		vec<-c(vec,sum(actual_val>mm[,i] & actual_val<mm[,i+1])/length(actual_val))
		i<-i+2
	}
	
	stats$q<-c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95) 
	stats$coverage<-rev(vec)
	
	return(stats)
}

model.coverage<-function(mod.pred,data.est,actual_est){
	library(boot)

	
	#estimation ###########################################################################
	index.est = inla.stack.index(stack,"est")$data
	stats<-list()
	
	qr<-function(q){return(cbind(.5*(1.-q),1-(.5*(1.-q))))}
	quants<-c(qr(0.95),qr(0.85),qr(0.75),qr(0.65),qr(0.55),qr(0.45),qr(0.35),qr(0.25),qr(0.15),qr(0.05))
	mm<-matrix(ncol=length(quants),nrow=length(index.est))
	for(i in index.est){
		mm[i,]<-inla.qmarginal(quants,mod.pred$marginals.linear.predictor[[i]])
	}
	mm<-inv.logit(mm)
	
	vec<-c()
	i<-1
	while(i < (length(quants))){
		vec<-c(vec,sum(actual_est>mm[,i] & actual_est<mm[,i+1])/length(actual_est))
		i<-i+2
	}
	
	stats$q<-c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95) 
	stats$coverage<-rev(vec)
	
	return(stats)
}


hpd.model.coverage<-function(mod.pred,data.est){
	library(boot)
	index.est = inla.stack.index(stack,"est")$data
	stats<-list()
	actual_est<-data.est[,'response']/data.est[,'examined']

	quants<-c((0.95),(0.85),(0.75),(0.65),(0.55),(0.45),(0.35),(0.25),(0.15),(0.05))

	mm_low<-matrix(ncol=length(quants),nrow=length(index.est))
	mm_high<-matrix(ncol=length(quants),nrow=length(index.est))

	for(i in index.est){
		mm<-inla.hpdmarginal(quants,mod.pred$marginals.linear.predictor[[i]])
		mm_low[i,]<-as.numeric(mm[,1])
		mm_high[i,]<-as.numeric(mm[,2])	
	}
	mm_low<-inv.logit(mm_low)
	mm_high<-inv.logit(mm_high)

	vec<-c()

	for(i in 1:length(quants)){
		vec<-c(vec,sum(actual_est>mm_low[,i] & actual_est<mm_high[,i])/length(actual_est))
	}
	
	stats$q<-c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95) 
	stats$coverage<-rev(vec)
	
	return(stats)
}

 
temporalInfoqtr<-function(f,dateEnd){
 library(zoo)
 start_dates<-as.yearmon(paste(f[, "year_start"], "-", f[, "month_start"],sep=""))
 end_dates<-as.yearmon(paste(f[, "year_end"], "-", f[, "month_end"],sep=""))
 mid_dates<-(start_dates+end_dates)/2
 fix_date<-as.yearmon(dateEnd)
 times<-as.yearqtr(fix_date)-as.yearqtr(mid_dates)
 return(times)
}  

temporalInfo<-function(f,dateEnd){
 library(zoo)
 start_dates<-as.yearmon(paste(f[, "year_start"], "-", f[, "month_start"],sep=""))
 end_dates<-as.yearmon(paste(f[, "year_end"], "-", f[, "month_end"],sep=""))
 mid_dates<-(start_dates+end_dates)/2
 fix_date<-as.yearmon(dateEnd)
 times<-fix_date-mid_dates
 return(times)
}  

splitdf <- function(dataframe, seed=NULL,Tp) {
	if (!is.null(seed)) set.seed(seed)	
	index <- 1:nrow(dataframe)	
	trainindex <- sample(index, round(Tp*nrow(dataframe)))	
	trainset <- dataframe[trainindex, ]	
	testset <- dataframe[-trainindex, ]	
	list(trainset=trainset,testset=testset)
}


xyz.to.ll <- function(cartesian)
{
	x<-cartesian[,1]
	y<-cartesian[,2]
	z<-cartesian[,3]

    r = sqrt(x*x + y*y+ z*z); 
    lat = asin(z/r)*(180/pi);
    lon = atan2(y, x)*(180/pi);
    return(cbind(lon,lat));
}



load.data<-function(loc,raster.fnames,fnames,reg,country,wh,isSpatial,qtr,reposition,dateEnd){
	dateEnd<-as.yearmon(dateEnd)
	f<-read.csv(loc)
	if(wh=='ALL'){ # use all the data 
	
	}else{
		if(reg & is.numeric(wh) ){	#use a region
			if('map_region' %in% colnames(f)){
				f<-f[f[,'mcmc_region_id']==wh,]
				print(paste('Using region',wh))
			} else {
				r=raster('/home/drive/mastergrids/cleandata/mcmcpfreg/mcmcpfreg5k_CLEAN.flt')
				which.region<-extract(r,cbind(f$longitude,f$latitude))	
				f<-f[which.region==wh,]		
				print(paste('Using region',wh))			
			}
		}
		if(country){
			cn<-raster('/home/likewise-open/ZOO/bras2280/Dengue/Rasters/cn5k.flt')
			gauls<-extract(cn,cbind(f[,'longitude'],f[,'latitude']))
			index<-gauls==wh
			f<-f[index,]
			print(paste('Using country',wh))
		}
		if(!is.numeric(wh) & wh!='ALL'){
			f<-f[f[,'map_region']==wh,]
			print(paste('Using region',wh))
		}
	}
	if(isSpatial){
		print(paste('CHOICE: spatial only data'))
	}else{
		print(paste('CHOICE: spatio-temporal data'))	
	}
	if(reposition){
	 	print(paste('CHOICE: center lat lons to 5k pixel'))			
	}else {
		print(paste('CHOICE: do not center lat lons to 5k pixel - use raw data'))	
	}
	
	
	
	# do not change the data
	if(qtr==TRUE){
		print(paste('CHOICE: aggregate monthly data to quaterly, backward referenced from',dateEnd))				
		months=temporalInfoqtr(f,dateEnd)	# get qtr
	} else {
		print(paste('CHOICE: using monthly data, backward referenced from',dateEnd))			
		months=temporalInfo(f,dateEnd)	# get months
	}
	f<-cbind(f,months)
	# for spatial only data snap points in 5km pixel to centroid and then aggregate number examined and positive
	r<-raster('/home/likewise-open/ZOO/bras2280/Sensitivity/Rasters/ls5k.flt')
	ce<-cellFromXY(r,f[,c('longitude','latitude')])  # get cell numbers
	cexy<-xyFromCell(r,ce)	# get xy from cells
	catData<-rep(NA,length(ce))		
	dup<-duplicated(catData)
	# concatenate lat-lon-months for unique entry searching
	for(i in 1:length(ce)){
		if(isSpatial){
			if(reposition){
				 catData[i]<-paste(cexy[i,1],cexy[i,2],sep=':')				
			 } else {
			 	 catData[i]<-paste(f[i,c('longitude')],f[i,c('latitude')],sep=':')				
			 }
		}else{
			if(reposition){
				 catData[i]<-paste(cexy[i,1],cexy[i,2],months[i],sep=':')
			} else {
			 	 catData[i]<-paste(f[i,c('longitude')],f[i,c('latitude')],months[i],sep=':')				
			 }				 
		}
	}

	un<-unique(catData)		#find points with unique lat,lon and time
	#### required columns of f
	f<-f[,c('longitude','latitude','positive','examined','lower_age','upper_age','months')]
	fcolapse<-matrix(nrow=length(un),ncol=ncol(f))	#create colapse matrix to sum aggregated data
	colnames(fcolapse)<-colnames(f)				#label			
	f[,'longitude']<-cexy[,1]	# add fixed lons
	f[,'latitude']<-cexy[,2]	#add fixed lats
	for(i in 1:length(un)){
		tmp<-f[catData%in%un[i],]
		fcolapse[i,'longitude']<-tmp[1,'longitude']
		fcolapse[i,'latitude']<-tmp[1,'latitude']		
		fcolapse[i,'months']<-tmp[1,'months']													
		fcolapse[i,'positive']<-as.integer(sum(tmp[,'positive']))
		fcolapse[i,'examined']<-as.integer(sum(tmp[,'examined']))				
		fcolapse[i,'lower_age']<-as.integer(mean(tmp[,'lower_age']))
		fcolapse[i,'upper_age']<-as.integer(mean(tmp[,'upper_age']))	

	}

	f<-fcolapse
	xyz<-ll.to.xyz(f[,c('longitude','latitude')])
	bindata<-f[,c('positive','examined')]
	cov<-extract(stack(raster.fnames),f[,c('longitude','latitude')])
	colnames(cov)<-fnames
	
	data<-list()
	data$x<-xyz[,'x']
	data$y<-xyz[,'y']
	data$z<-xyz[,'z']
	data$response<-bindata[,'positive']
	data$examined<-bindata[,'examined']
	data$lower_age=f[,"lower_age"]
	data$upper_age=f[,"upper_age"]
	data$months=f[,'months']
	outdata<-cbind(as.data.frame(data),as.data.frame(cov))
	print(paste('This dataset contains - -',nrow(outdata),'points and',ncol(cov),'covariates'))
	return(outdata)
}

load.evi<-function(data){
		library(raster)
		library(zoo)
		library(mgcv)
		library(cvTools)
		f<-read.csv('~/Desktop/Space Time Malaria/PR data May 2013.csv')
		f<-f[f$mcmc_region_id==1,]
		#f<-f[f$map_region=='Africa+',]
		
		 start_dates<-as.yearmon(paste(f[, "year_start"], "-", f[, "month_start"],sep=""))
		 end_dates<-as.yearmon(paste(f[, "year_end"], "-", f[, "month_end"],sep=""))
		 mid_dates<-(start_dates+end_dates)/2
		 
		 cutoffdata<-as.yearmon('2000-04')
		 f<-f[mid_dates>cutoffdata,]
		 mid_dates<-mid_dates[mid_dates>cutoffdata]
		 
		 split<-strsplit(as.character(mid_dates),split=" ")
		 years<-rep(NA,length(mid_dates))
		 months<-rep(NA,length(mid_dates))
		 
		 for(i in 1:length(mid_dates)){
			tmp<-split[[i]]
			months[i]<-tmp[1]
			years[i]<-tmp[2]
		 
		 }
		 
		 charmon<- sort(unique(months))
		 numericmon<-c('04','08','12','02','01','07','06','03','05','11','10','09')
		 numericyear<-as.character(2000:2012)
		 for(i in 1:length(charmon)){
			months[months==charmon[i]]<-as.character(numericmon[i])
		 }
		 response<-log((as.numeric(f$positive)+0.5)/(as.numeric(f$examined)-as.numeric(f$positive)+0.5))
		 d<-data.frame(f$longitude,f$latitude,f$positive,f$examined,months,years,1:length(months),response=response)
		colnames(d)<-c("longitude","latitude","positive","examined","months","years",'id','response')
		
		evi<-rep(NA,nrow(d))
		response<-rep(NA,nrow(d))
		i=1
		sortedmon<-sort(numericmon)
		for(y in numericyear){
			j=1
			for(m in sortedmon){
				if(paste(m,y,sep='-')!='01-2000'){
					tmp=d[d$months==m & d$years==y,]
					if(nrow(tmp)!=0){
						if(m=='01'){
						   lm='12'		
						   ly=as.character(as.numeric(y)-1)
						}else {
						  lm=sortedmon[j-1]			
						  ly=y
						}
						r<-raster(paste('/Volumes/map_data/users/sam/EVI_5km_data/',ly,'.',lm,'.EVI.tif',sep=""))
						evi[tmp$id]<-extract(r,cbind(as.numeric(tmp$longitude),as.numeric(tmp$latitude)))
						response[tmp$id]<-tmp$response
						print(paste(m,y,'has been lagged by 1 to',lm,ly))
					}
				}
				j=j+1
			}
			i=i+1
		}
		evi.store1=evi
		response.store1=response


}


pred.data<-function(raster.fnames,r,reg){
	s<-stack(raster.fnames)
	r[r!=reg]=NA
	r<-trim(r) 
	pred_val<-getValues(r)
	w<-is.na(pred_val)
	pred_locs<-xyFromCell(r,1:ncell(r))  
	pred_locs<-pred_locs[!w,] 
	pred_data<-extract(s,pred_locs)
	colnames(pred_locs)<-c('longitude','latitude')
	pred_locs<-ll.to.xyz(pred_locs)


}

gc.dist<-function(r){
	return(6356.7523142*2*asin(r/2))
}

distances<-function(xyz){
	distances<-rep(0,nrow(xyz))
	for(i in 1:nrow(xyz)){
		x2<-(xyz[i,'x']-xyz[-i,'x'])^2
		y2<-(xyz[i,'y']-xyz[-i,'y'])^2
		z2<-(xyz[i,'z']-xyz[-i,'z'])^2
		distances[i]<-mean(sqrt(x2+y2+z2))
	
	}
	norm.dist<-( distances-min(distances))/(max(distances)-min(distances))
	return(norm.dist)
}


distances.matrix<-function(xyz){
	distances<-matrix(nrow=nrow(xyz),ncol=nrow(xyz))
	for(i in 1:nrow(xyz)){
		x2<-(xyz[i,1]-xyz[,1])^2
		y2<-(xyz[i,2]-xyz[,2])^2
		z2<-(xyz[i,3]-xyz[,3])^2
		distances[i,]<-sqrt(x2+y2+z2)
	
	}
	return(distances)
}


create.est.stack.fit<-function(data.est,mesh,knots,interval){
	mesh1d=inla.mesh.1d(knots,interval=interval,degree=2, boundary='free')
	
	
	est.cov<-create.list.cov(data.est)
	
	A.est =
	  inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'months'],group.mesh=mesh1d)
	  
		#-- Create index matrix --#
	field.indices =
	  inla.spde.make.index("field", n.spde=mesh$n,n.group=mesh1d$m)
	
	 
	iid.indices.est =
	  inla.spde.make.index("iid", n.spde=nrow(data.est))
	  
	  stack.est =
	inla.stack(data=list(response=data.est$response,examined=data.est$examined),
			   A=list(A.est,1),
			   effects=
				 list(c(field.indices,
						list(Intercept=1)),
						c(est.cov,iid.indices.est)
					  ),
			   tag="est", remove.unused=TRUE)
		return(stack.est)
}

create.est.stack.pred<-function(data.est,mesh,knots,interval){
		mesh1d=inla.mesh.1d(knots,interval=interval,degree=2, boundary='free')
		
		l=mesh$n*mesh1d$m
		mat=matrix(data=NA,nrow=l,ncol=ncol(data.est))
		colnames(mat)<-colnames(data.est)
		mat<-as.data.frame(mat)
		data.est<-as.data.frame(data.est)
		mat[1:nrow(data.est),]<-data.est
		data.est<-mat
		num.dates.missing<-sum(is.na(data.est$months))
		data.est$months[is.na(data.est$months)]=sample(loc,num.dates.missing,replace=TRUE)
		
		est.cov<-create.list.cov(data.est)
		A.est =
		  inla.spde.make.A(mesh, loc=cbind(data.est[,'x'],data.est[,'y'],data.est[,'z']),group=data.est[,'months'],group.mesh=mesh1d)
		  		
			#-- Create index matrix --#
		field.indices =
		  inla.spde.make.index("field", n.spde=mesh$n,n.group=mesh1d$m)
		
		
		iid.indices.est=as.data.frame(1:length(field.indices$field))
		colnames(iid.indices.est)<-'iid'
			
		
		  stack.est =
		inla.stack(data=list(response=data.est$response,examined=data.est$examined),
				   A=list(A.est,1),
				   effects=
					 list(c(field.indices,
							list(Intercept=1)),
							c(est.cov,iid.indices.est)
						  ),
				   tag="est",remove.unused=FALSE)
			return(stack.est)	   
}

create.val.stack.fit<-function(data.val,mesh,knots,interval){
	mesh1d=inla.mesh.1d(knots,interval=interval,degree=2, boundary='free')
	
	val.cov<-create.list.cov(data.val)
	
	A.val =
	  inla.spde.make.A(mesh, loc=cbind(data.val[,'x'],data.val[,'y'],data.val[,'z']),group=data.val[,'months'],group.mesh=mesh1d)
	  
		#-- Create index matrix --#
	field.indices =
	  inla.spde.make.index("field", n.spde=mesh$n,n.group=mesh1d$m)
	
	 
	iid.indices.val =
	  inla.spde.make.index("iid", n.spde=nrow(data.val))
	  
	  stack.val =
	inla.stack(data=list(response=NA,examined=NA),
			   A=list(A.val,1),
			   effects=
				 list(c(field.indices,
						list(Intercept=1)),
						c(val.cov,iid.indices.val)
					  ),
			   tag="val", remove.unused=TRUE)
			   
		return(stack.val)
}