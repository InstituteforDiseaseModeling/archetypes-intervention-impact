## functions from "Create_database_ITN... etc"
library(zoo)
library(raster)
library(VGAM)

# functions: 
# -emplogit2: 2-input empirical logit
# -reposition.points: adjust lat-longs (for what reason?)
# -aggregate.data: sum to pixel level
# -calc_access_matrix and get.access.mat: access calculations from (?)



######################################### FUNCTIONS ###########################################################

emplogit2<-function(Y,N) {top=Y+0.5;bottom=N-Y+0.5;return(log(top/bottom))} # empirical logit two numbers

reposition.points<-function(ls,points,square){
  
  rm<-extract(ls,cbind(points[,'lon'],points[,'lat']))	# get all point values
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
  
  # accdate1: time at which to calculate access
  # func0: list of functions defining the probability a house of size n has no nets
  # func1: list of mean # of nets per household for a house of size n
  # hh: vector of proportion of houses at different sizes (1-10+)
  # nc: number of columns (todo: calculate from length of func0, func1, or hh)
  # nr: number of rows?
  
  # returns: a vector of access metrics.
  
  # find probability of no net and mean # of nets across household size for the time in question
  k1<-k0<-rep(NA,10)
  for(j in 1:length(func0)){
    k1[j]<-func1[[j]](accdate1)
    k0[j]<-func0[[j]](accdate1)
  }
  
  # fill na's (when would k0 be null?)
  k0[is.nan(k0)]<-0
  
  # make a 41x10 matrix
  matModel<-matrix(data=0,nrow=nr+1,ncol=nc)
  cnames<-colnames(matModel)<-1:nc
  rnames<-rownames(matModel)<-0:nr
  
  # first row of matrix: household size distribution * probability of no net (weighted probability of no net ownership)
  # don't think these values are ever used?
  
  matModel[1,]<-hh*k0
  
  # "remaining": household size distribution * probability of having a net (weighted probability of net ownership)
  remaining<-hh*(1-k0)
  for(i in 1:ncol(matModel)){
    l<-k1[i] # mean # of nets in household of size i
    
    # Find the probability of a household of size i having 1-nr nets, assuming a zero-truncated poisson distribution with mean equal 
    # to the stock and flow mean. Multiply this by the (weighted) probability of having a net.
    # Each column of matModel (not including the first row) sums to the corresponding value of "remaining".
    matModel[2:nrow(matModel),i]<-dpospois(1:(nrow(matModel)-1),l)*remaining[i]
  }	
  
  cnames<-colnames(matModel)<-1:nc
  rnames<-rownames(matModel)<-0:nr
  
  # calculate access
  ind3mat<-matrix(data=0,nrow=nrow(matModel),ncol=ncol(matModel))
  ind3tot<-matrix(data=0,nrow=nrow(matModel),ncol=ncol(matModel))
  
  for(x in 1:ncol(matModel)){ # x is household size
    for(y in 1:nrow(matModel)){ # y is net count
      ind3mat[y,x]<-2*(y-1)*matModel[y,x] # prop of people with access (defined as 1-between-2) to y nets
      ind3tot[y,x]<-x*matModel[y,x] # prop of people in this net-hhsize combo 
      if(ind3mat[y,x]>ind3tot[y,x]) ind3mat[y,x]=ind3tot[y,x] #cap on estimate 
    }
  }	
  
  acc<-rep(NA,10)
  
  # access by household
  for(i in 1:10){
    acc[i]<-sum(ind3mat[,i])/sum(ind3tot[,i])
  }
  
  # overall access
  acc<-c(acc,sum(ind3mat)/sum(ind3tot))
  return(acc)
}

