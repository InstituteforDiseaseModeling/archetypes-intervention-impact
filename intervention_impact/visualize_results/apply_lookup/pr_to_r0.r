rb = function(N=1){rbeta(N,82,68)}

rc = function(N=1){rnorm(N,0.161,0.0092)} 

rD = function(N=1){rnorm(N,33.8,1.94)} 

rkappa = function(N=1){rgamma(N, shape=9.36, scale=0.0077)}

rS = function(N=1){rlnorm(N, 0.8, 0.765)}

ralpha = function(N=1){ rnorm(N, 4, 1) }

rPR2EIR = function(PR, N=1){
  mu = -1.768 + 7.247*PR; sigma = 1.281
  rlnorm(N, mu, sigma)
}

rPR2R.fac = function(N=1){
  kappa = rkappa(N) 
  rb(N)*rD(N)*(1+ralpha(N))*kappa*(1+rS(N)*kappa)/kappa/365
}

rPR2R.F1 = function(PR, N=1){
  S = rS(N) 
  #alpha = ralpha(PR,N) 
  alpha = ralpha(N) 
  c = rc(N)
  ((1-PR)^(-alpha) -1) / (1-(1-PR)^(1+alpha)) *
    (1+alpha)/alpha*(1+S*c*(1-(1-PR)^(1+alpha)))
} 

rPR2R.F2 = function(PR, N=1){rPR2EIR(PR,N)*rPR2R.fac(N)} 

weight = function(PR){exp(-80*PR^4)} 

rPR2R=function(PR, N=1){
  W = weight(PR)  
  rPR2R.F1(PR, N)*W + rPR2R.F2(PR,N)*(1-W)  
}

PR2R = function(PR, N=10000){
  median(rPR2R(PR,N),na.rm=T)
}

R2spline <- function(){
  PR<-seq(0,1,length=1002)[2:1001]
  Ro<-rep(NA,length(PR))
  Ro <- unlist(lapply(PR, PR2R))
  splObj<-smooth.spline(PR, Ro)
  plot(PR,Ro)
  Xpred<-seq(0,1,length=100)
  Ypred<-as.vector(as.matrix(predict(splObj,data.frame("X"=Xpred))$y))
  lines(Xpred,Ypred,col=2)
  return(splObj)
}


