###############################################################################################################
## 03_05_general_functions.r
## Amelia Bertozzi-Villa
## May 2019
## 
## Functions to accompany 03_access_dev.r 

## -emplogit: 1-value empirial logit. 
## -other_emplogit: compare to emplogit
## -emplogit2: 2-input empirical logit
## - ll.to.xyz: convert a set of lat-longs to cartesian gridpoints.

## -reposition.points: adjust impossibly-placed lat-longs
## -aggregate.data: sum to pixel level (TODO: fix column bug)
## -calc_access_matrix: calculate access metrics by household size from stock and flow outputs
##############################################################################################################


emplogit<-function(Y,tol){
  # Y: value to transform
  # tol: tolerance value to prevent zeros
  top=Y*tol+0.5
  bottom=tol*(1-Y)+0.5
  return(log(top/bottom))
}

other_emplogit <- function(Y, tol){
  top = tol + Y
  bottom = 1 - Y + tol
  return(log(top/bottom))
}

emplogit2<-function(Y,N){
  # approximation of a log odds
  # Y: # of occurrences of interest
  # N: # of tries
  top=Y+0.5
  bottom=N-Y+0.5
  return(log(top/bottom))
}


# Inverse Hyperbolic sin transform
ihs <- function(x, theta){  # function to IHS transform
  return(asinh(theta * x)/theta) 
}

# Inverse of the inverse hyperbolic sin transform
inv_ihs <- function(x, theta){
  (1/theta)*sinh(theta * x)
}

# Inverse hyperbolic sin transform-- log-likelihood
ihs_loglik <- function(theta,x){
  
  n <- length(x)
  xt <- ihs(x, theta)
  
  log.lik <- -n*log(sum((xt - mean(xt))^2))- sum(log(1+theta^2*x^2))
  return(log.lik)
}

ll_to_xyz<-function(ll){
  
  ## ll: data.table with columns "row_id", "longitude", "latitude"
  ll <- ll[, list(row_id, longitude, latitude,
                  longitude_rad=longitude*(pi/180),
                  latitude_rad=latitude*(pi/180))]
  
  xyz <- ll[, list(row_id,
                   x=cos(latitude_rad) * cos(longitude_rad),
                   y=cos(latitude_rad) * sin(longitude_rad),
                   z=sin(latitude_rad))]
  
  return(xyz)
}