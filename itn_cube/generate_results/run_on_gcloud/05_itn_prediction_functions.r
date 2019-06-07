


package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "data.table", "rgdal", "INLA", "RColorBrewer", "VGAM", "maptools", "rgeos", "png", "sp", "doParallel", "SDMTools"))



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


get.pred.locs<-function(in_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/input_data_archive/general"){
  cn<-raster(file.path(in_dir, 'african_cn5km_2013_no_disputes.tif')) #load raster
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



get.pred.covariates.dynamic<-function(year,month, cov_dir="/mnt/data/input/gs/map_data_z/cubes_5km", 
                                      joint_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/input_data_archive/general"){
  
  dynamic_covnames <- c("LST_day/mean", "LST_night/mean", "EVI/mean", "TCW/mean", "TSI/mean")
  #' foldersd<-c('/home/drive/cubes/5km/LST_day/mean/',
  #'             #'/home/drive/cubes/5km/LST_delta/mean/',
  #'             '/home/drive/cubes/5km/LST_night/mean/',
  #'             '/home/drive/cubes/5km/EVI/mean/',
  #'             #'/home/drive/cubes/5km/TCB/mean/',
  #'             '/home/drive/cubes/5km/TCW/mean/',
  #'             '/home/drive/cubes/5km/TSI/mean/')
  
  l=length(dynamic_covnames)
  
  cn<-raster(file.path(joint_dir, 'african_cn5km_2013_no_disputes.tif')) #load raster
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
  
  ncores <- detectCores()
  registerDoParallel(ncores-2)
  covs.list.dyn<-foreach(i=1:l,.combine='cbind') %dopar% { # loop through unique names
    r=raster(file.path(cov_dir, dynamic_covnames[i], paste0(year,'.',month,'.mean.tif')))
    NAvalue(r)=-9999
    return(r[index])
  }
  colnames(covs.list.dyn)<-c('lst_day','lst_night','evy','tcw','tsi')
  
  return(covs.list.dyn)
}

get.pred.covariates.yearonly<-function(year,cov_dir="/mnt/data/input/gs/map_data_z/cubes_5km", 
                                       joint_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/input_data_archive/general/"){
  # foldery<-c('/home/drive/cubes/5km/IGBP_Landcover/Fraction/','/home/drive/cubes/5km/AfriPop/')
  cn<-raster(file.path(joint_dir, 'african_cn5km_2013_no_disputes.tif')) #load raster
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
  
  ncores <- detectCores()
  registerDoParallel(ncores-2)
  covs.list.year<-foreach(i=1:l,.combine='cbind') %dopar% { # loop through unique names
    # no land cover for 2013
    if(year>2012){
      r=raster(file.path(cov_dir, "IGBP_Landcover/Fraction", paste0('2012','.fraction.class.',i-1,'.tif'))) 
    }else{
      r=raster(file.path(cov_dir, "IGBP_Landcover/Fraction", paste0(year,'.fraction.class.',i-1,'.tif')))
    }		
    NAvalue(r)=-9999
    return(r[index])
  }
  covs.list.year<-covs.list.year[,-14]	#  # remove landcover 13  - Urban and built-up   - for collinearity with population
  covnames<-0:16
  covnames<-covnames[-14]	
  r=raster(file.path(cov_dir, "AfriPop", paste0(year,'.total.population.tif')))
  NAvalue(r)=-9999
  pop<-r[index]
  covs.list.year<-cbind(covs.list.year,pop)
  colnames(covs.list.year)<-c(paste0('landcover',covnames),'populatopn')
  
  return(covs.list.year)	
}

get.pred.covariates.static<-function(cov_dir="/mnt/data/input/gs/map_data_z/cubes_5km", 
                                     joint_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/input_data_archive/general"){
  
  #### static covariates
  #### folders
  #' folderss<-c('/home/drive/cubes/5km/Topographic/Africa_TMI_90m.mean.tif',
  #'             '/home/drive/cubes/5km/Topographic/Africa_SRTM_90m.mean.tif',
  #'             '/home/drive/cubes/5km/Topographic/Africa_slope_90m.mean.tif',
  #'             #	'/home/drive/cubes/5km/Topographic/Africa_FA_90m.mean.tif',
  #'             '/home/drive/cubes/5km/Seasonality/pf_seasonality.tif',
  #'             '/home/drive/cubes/5km/Poverty/PET_1950-2000_5km.mean.tif',
  #'             '/home/drive/cubes/5km/Poverty/AI_1950-2000_5km.mean.tif',
  #'             '/home/drive/cubes/5km/Poverty/accessibility_50k_5km.mean.tif',
  #'             #'/home/drive/cubes/5km/worldclim/prec57a0.tif',
  #'             '/home/drive/cubes/5km/Poverty/viirs_nighttime_5km.mean.tif')
  #' #'/home/drive/cubes/5km/Poverty/DMSP_F18_5km.mean.tif')
  
  static_fnames <- c("Topographic/Africa_TMI_90m.mean.tif",
                     "Topographic/Africa_SRTM_90m.mean.tif",
                     "Topographic/Africa_slope_90m.mean.tif",
                     "Seasonality/pf_seasonality.tif",
                     "Poverty/PET_1950-2000_5km.mean.tif",
                     "Poverty/AI_1950-2000_5km.mean.tif",
                     "Poverty/accessibility_50k_5km.mean.tif",
                     "Poverty/viirs_nighttime_5km.mean.tif")
  
  static_covnames <- file.path(cov_dir, static_fnames)
  
  cn<-raster(file.path(joint_dir, 'african_cn5km_2013_no_disputes.tif')) #load raster
  NAvalue(cn)=-9999
  pred_val<-getValues(cn)#get values again
  w<-is.na(pred_val) #find NAs again
  index<-1:length(w) 
  index<-index[!w]
  pred_locs<-xyFromCell(cn,1:ncell(cn))  #get prediction locations
  pred_locs<-pred_locs[!w,] #remove NA cells
  
  st<-stack(static_covnames)
  NAvalue(st)=-9999
  
  static.covs<-st[index]
}


# from later: plotting 
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

rel_gain_hist<-function(in_dir){ # previously: in_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/predictions"
  h<-c()
  P=raster::stack(paste0(in_dir, '/ITN_',2000:2016,'.ACC.tif'))
  P5=stack(paste0(in_dir, '/ITN_',2000:2016,'.USE.tif'))
  ranges<-matrix(nrow=nlayers(P5),ncol=2)
  for(i in 1:nlayers(P5)){
    tmp<-relative_gain(P5[[i]],P[[i]])
    tmp<-getValues(tmp)
    tmp<-tmp[!is.na(tmp)]
    h<-c(h,tmp)
  }
  return(h)
  
}

rel_gain_range<-function(in_dir){ # previously: in_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/predictions"
  
  P=raster::stack(paste0(in_dir, '/ITN_',2000:2016,'.ACC.tif'))
  
  P5=stack(paste0(in_dir, '/ITN_',2000:2016,'.USE.tif'))
  ranges<-matrix(nrow=nlayers(P5),ncol=2)
  for(i in 1:nlayers(P5)){
    tmp<-relative_gain(P5[[i]],P[[i]])
    tmp<-range(getValues(tmp),na.rm=T)
    ranges[i,]<-tmp
  }
  return(c(min(ranges[,1]),max(ranges[,2])))
  
}

dev_range<-function(in_dir){ # previously: in_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/predictions"
  P=raster::stack(paste0(in_dir, '/ITN_',2000:2016,'.DEV.tif'))
  ranges<-matrix(nrow=nlayers(P),ncol=2)
  for(i in 1:nlayers(P)){
    ranges[i,]<-ranges[i,]<-cbind(cellStats(P[[i]],stat='min'),cellStats(P[[i]],stat='max'))
  }
  return(c(min(ranges[,1]),max(ranges[,2])))
  
}

gap_range<-function(in_dir){ # previously: in_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/predictions"
  P=raster::stack(paste0(in_dir, '/ITN_',2000:2016,'.GAP.tif'))
  ranges<-matrix(nrow=nlayers(P),ncol=2)
  for(i in 1:nlayers(P)){
    ranges[i,]<-ranges[i,]<-cbind(cellStats(P[[i]],stat='min'),cellStats(P[[i]],stat='max'))
  }
  return(c(min(ranges[,1]),max(ranges[,2])))
  
}


## what is cn? country? 
get.cn.estimates<-function(r,pop, in_dir="/mnt/data/input/gs/map_data_z/users/amelia/itn_cube/input_data_archive"){
  database<-read.csv(file.path(in_dir, 'for_plots_and_tables/country_table.csv'))
  cn<-raster(file.path(in_dir, 'general/african_cn5km_2013_no_disputes.tif')) #load raster
  NAvalue(cn)=-9999
  try(if(extent(cn)!=extent(r)) stop("extents do not match"))
  limits<-raster(file.path(in_dir, 'for_plots_and_tables/Pf_Limits/Pf_limits.tif'))
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
