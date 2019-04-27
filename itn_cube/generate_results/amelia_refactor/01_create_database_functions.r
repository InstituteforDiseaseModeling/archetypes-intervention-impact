###############################################################################################################
## 01_create_database_functions.r
## Amelia Bertozzi-Villa
## April 2019
## 
## Functions to accompany create_database_refactored.r-- prepping survey data to go into the ITN cube model. 
## -emplogit2: 2-input empirical logit
## -reposition.points: adjust impossibly-placed lat-longs
## -aggregate.data: sum to pixel level (TODO: fix column bug)
## -calc_access_matrix: calculate access metrics by household size from stock and flow outputs
##############################################################################################################

emplogit2<-function(Y,N){
  # approximation of a log odds
  # Y: # of occurrences of interest
  # N: # of tries
  top=Y+0.5
  bottom=N-Y+0.5
  return(log(top/bottom))
  }

reposition.points<-function(nat_raster, points, radius=4){
  
  # nat_raster: raster layer of country boundaries
  # points: data.table including the columns "lat" and "long" for each cluster location
  # radius: cell count "radius" around which to search for a point to reposition
  
  matrix_side_len <- radius*2 + 1
  rm<-extract(nat_raster, points[, list(lon,lat)])	# get all point values
  index<-which(is.na(rm))
  
  # loop through the points that need to be repositioned
  for(x in index){
    
    st<-paste0('<', 'Survey-year: ', points$Survey[x],'-',points$year[x],'>')
    error <- paste0("--> ", st, " Point ",  points[x,'lat'], points[x,'lon']," is not positioned properly")
    
    # find cell 4 rows up and 4 columns left of the point of interest (for matrix contstruction)
    row1 =  rowFromY(nat_raster, points[[x,'lat']])-radius 
    col1 = colFromX(nat_raster, points[[x,'lon']])-radius 
    
    # if these are valid row values, attempt to reposition
    if(!is.na(col1) & !is.na(row1)){ 
      
      # construct a block of values with a radius of "radius" around the point to be repositioned
      mat<-raster::getValuesBlock(nat_raster, row=row1, nrows=matrix_side_len, col=col1, ncols=matrix_side_len) 

      # if any of the points *near* the original points exist, then reposition point.
      if(any(!is.na(mat))){	
        
        mat <- matrix(mat, nrow=matrix_side_len, ncol=matrix_side_len, byrow=TRUE)
        matcells<-mat
        
        # get cell counts for each point in this
        for(i in 1:dim(mat)[1]){
          for(j in 1:dim(mat)[2]){
            matcells[i,j]=as.numeric(cellFromRowCol(nat_raster,(row1-1+i),(col1-1+j)))	
          }
        }
        
        # get raw lat lons for each cell above
        matcoords<-xyFromCell(nat_raster,matcells) 
        matcells<-as.vector(matcells)
        
        # keep only non-null nearby pixels
        matcells<-matcells[!is.na(as.vector(mat))]
        matcoords<-matcoords[!is.na(as.vector(mat)),]
        
        # find distances between current point and non-null points 
        dist<-pointDistance(cbind(points[x,'lon'],points[x,'lat']),matcoords,longlat=T)	# find the closest non na valid cell
        
        # assign point a new lat-long equal to the closest non-null point
        newp<-xyFromCell(nat_raster,matcells[which.min(dist)])
        points[x,'lon']<-newp[,'x']
        points[x,'lat']<-newp[,'y']
        
        print(paste('--> Point', st ,'has been repositioned at distance of', min(dist)/1000, "km", "at(",points[x,'lon'],",",points[x,'lat'],")"))  	
        
      } else {
        print(error)
      }
    } else{
      print(error)
      } 
  }
  return(points)
  
}


aggregate.data<-function(nat_raster, data){
  
  # snap points in each 5km pixel to centroid and then aggregate number with access and number in household
  # TODO: why only aggregate those two numbers?
  # nat_raster: raster layer of country boundaries
  # data: data.table to aggregate

  ce<-cellFromXY(nat_raster, data[, list(lon, lat)])  # get cell numbers
  data$cellnumber<-ce # update the cell numbers
  cexy<-xyFromCell(nat_raster, ce)	# get xy from cells
  # update lat/long values
  data[, lon:=cexy[,1]]
  data[, lat:=cexy[,2]]
  
  # within each cell: update P and N to the sum of clusters within the cell
  data[, P:= sum(P), by=list(cellnumber, yearqtr)] # sum those with access to a net
  data[, N:= sum(N), by=list(cellnumber, yearqtr)] # sum those who slept in house the night before 
  
  # Keep only the first cluster in each cell. 
  # TODO: update this so that instead all columns, not just P and N, get updated before dropping other rows.
  data[, groupid:= seq_len(.N), by=list(cellnumber, yearqtr)]
  data <- data[groupid==1]
  data[, groupid:=NULL]
  
  return(data)

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

