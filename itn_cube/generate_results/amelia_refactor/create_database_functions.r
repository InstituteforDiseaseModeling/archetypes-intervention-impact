## functions from "Create_database_ITN... etc"
library(zoo)
library(raster)
library(VGAM)
library(data.table)

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


calc_access_matrix<-function(access_date, func0, func1, hh_props, max_nets=40){
  
  # access_date: time at which to calculate access
  # func0: list of functions defining the probability a house of size n has no nets in this country
  # func1: list of functions defining mean # of nets per household for a house of size n in this country
  # hh_props: vector of proportion of houses at different sizes (1-10+) based on this survey
  # max_nets: this function will calulate the probability of household-level ownership for a net count up to this variable
  
  hh_bins <- length(hh_props)
  
  # returns: a vector of length hh_bins+1, where the first hh_bins entries refer to access for a given household size 
  # (1:hh_bins) and the last entry refers to survey-wise access
  
  # find probability of no net and mean # of nets across household size for the time in question

  
  net_stats <- data.table(hh_size=1:hh_bins,
                          prob_no_nets=sapply(1:hh_bins, function(hh_size){func0[[hh_size]](access_date)}),
                          mean_nets_per_hh=sapply(1:hh_bins, function(hh_size){func1[[hh_size]](access_date)}),
                          hh_prop=hh_props)
  net_stats[, weighted_prob_no_nets:=hh_prop*prob_no_nets]
  net_stats[, weighted_prob_any_net:=hh_prop*(1-prob_no_nets)]
  
  # Find the probability of a household of size i having 1:max_nets nets, 
  # assuming a zero-truncated poisson distribution with mean equal to the stock and flow mean. 
  # Multiply this by the weighted probability of having a net.
  
  net_dist <- lapply(net_stats$hh_size, function(this_hh_size){
    
    these_net_stats <- net_stats[hh_size==this_hh_size]
    
    no_hh_subset <- data.table(hh_size=this_hh_size,
                               net_count=0,
                               weighted_net_prob=these_net_stats$weighted_prob_no_nets)
    
    net_hh_subset <- data.table(hh_size=this_hh_size,
                                net_count=1:max_nets,
                                weighted_net_prob=dpospois(1:max_nets, 
                                                           these_net_stats$mean_nets_per_hh) * these_net_stats$weighted_prob_any_net
                                )
    
    return(rbind(no_hh_subset, net_hh_subset))
    
  })
  
  net_dist <- rbindlist(net_dist)
  
  # calculate access 
  net_dist[, prop_with_access:=net_count*2*weighted_net_prob]
  net_dist[, prop_in_bin:=hh_size*weighted_net_prob]
  net_dist[, prop_with_access:=pmin(prop_with_access, prop_in_bin)] # cap access proportion at proportion of population
  
  access_stats <- net_dist[, list(access=sum(prop_with_access)/sum(prop_in_bin)), by=hh_size] 
  all_access_stats <- sum(net_dist$prop_with_access)/sum(net_dist$prop_in_bin)
  
  if (sum(net_dist$weighted_net_prob)!=1){
    warning("Net probabilities do not sum properly!")
  }
  
  return(c(access_stats$access, all_access_stats))
  
}

