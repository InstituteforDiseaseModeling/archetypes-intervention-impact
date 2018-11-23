
## Amelia's first pass at parsing "Create_database_ITN access....r"

## sam: takes Harry's data, does aggregating, find national level means at quartery times 

library(zoo)
library(raster)
library(VGAM)
library(doParallel)

source("create_database_functions.r")

# Data loading, household-level access/use stats  ------------------------------------------------------------

# p0 & p1-- from stock & flow, nat'l time series of p0=p(hh has >0 nets) and p1=avg # of nets
load('/home/backup/ITNcube/prop0prop1.rData')
KEY=read.csv('/home/backup/ITNcube/KEY_080817.csv')

# Cout must be an object in prop0prop1.rData
Country.list<-Cout

# load household data, keep only those in country list
HH<-read.csv('/home/backup/ITNcube/ALL_HH_Data_20112017.csv',stringsAsFactors=F)
HH<-HH[HH$ISO3%in%Country.list,]

# Remove households with size zero or NA
zerohh<-is.na(HH$n.individuals.that.slept.in.surveyed.hhs)
HH<-HH[!zerohh,] 

zerohh<-(HH$n.individuals.that.slept.in.surveyed.hhs)==0
HH<-HH[!zerohh,] 

# update survey vector
Surveys<-as.character(unique(HH$Survey.hh)) 
Surveys<-Surveys[!is.na(Surveys)]

# find access (# with a net available) and use (# sleeping under net) per household 
times<-seq(0.0,17,0.25)+2000 # 2000 to 2017, in quarter-year intervals
a=HH$n.individuals.that.slept.in.surveyed.hhs # number in household
b=HH$n.ITN.per.hh # number of nets
d=b*2 # number covered by nets
d[d>a]=a[d>a] # bounding to not exceed the number in the house
e=HH$n.individuals.that.slept.under.ITN # use count

# Main loop: calculating access/gap for each household cluster  ------------------------------------------------------------
registerDoParallel(62)
output<-foreach(i=1:length(Surveys),.combine=rbind) %dopar% { #survey loop
  svy<-Surveys[i] # store survey
  
  cc<-HH[HH$Survey.hh==svy,'ISO3'] # get country
  p<-out[[which(Country.list==cc[1])]]
  p0<-p[,,1] # get p0 parameters (what are these?--abv)
  p1<-p[,,2] # get p1 parameters	
  
  func0<-list()
  func1<-list()
  # make these into continuous functions (via linear interpolation--abv)
  for(f in 1:ncol(p0)){
    eval(parse(text=paste0('func0[[',f,']]<-approxfun(times,p0[,',f,'])')))
  }
  for(f in 1:ncol(p0)){ # should this be p1? may not matter if they're the same length
    eval(parse(text=paste0('func1[[',f,']]<-approxfun(times,p1[,',f,'])')))
  }
  
  tmp_hh=HH[HH$Survey.hh==svy,] # subset hosehold 
  un<-unique(tmp_hh$Cluster.hh) # get unique cluster IDs for that household survey
  
  ### house hold size as measured by the number of people sleeping under a net the night before
  ## (not sure what this is doing-- abv)
  hh_val<-rep(NA,100) 
  tmp_hh$sample.w<-tmp_hh$sample.w/sum(tmp_hh$sample.w) # get sample weights
  for(xx in 1:100){
    hh_val[xx]<-sum(tmp_hh[tmp_hh$n.individuals.that.slept.in.surveyed.hhs==xx,'sample.w'])	# proportion in a given household size	
  }
  hh_size<-hh_val[1:10] # 1:10
  hh_size[10]<-sum(hh_val[10:length(hh_val)])	# >10
  hh<-as.vector(as.numeric(hh_size))
  names(hh)<-1:10	
  
  #store matrix (1 row per cluster id, columns as described below)
  # p,n,t: data, means: means from stcok and flow
  fields<-c("Survey" ,"lat" , "lon",'year','P','N','Pu','T','Amean','Tmean','gap','gap2','gap3')
  tmp<-matrix(NA,nrow=length(un),ncol=length(fields))
  colnames(tmp)<-fields
  tmp<-as.data.frame(tmp)
  
  for(j in 1:length(un)){ #for each cluster
    wh<-HH$Survey.hh==svy & HH$Cluster.hh==un[j] # get the right HH values (what is wh? --abv)
    
    tmp[j,'Survey']<-svy
    tmp[j,'lat']<-mean(HH[wh,'latitude'],na.rm=TRUE)
    tmp[j,'lon']<-mean(HH[wh,'longitude'],na.rm=TRUE)
    
    tmp[j,'P']<-sum(d[wh]) #access cluster count
    tmp[j,'N']<-sum(a[wh]) #number in cluster count
    tmp[j,'year']<-mean(HH$year[wh])
    tmp[j,'Pu']<-sum(e[wh]) #number in cluster count (for USE--abv)
    tmp[j,'T']<-sum(b[wh]) #number in cluster count (number of NETS--abv)
    tmp[j,'gap']<-mean(((tmp[j,'P']/tmp[j,'N'])-(tmp[j,'Pu']/tmp[j,'N']))/(tmp[j,'P']/tmp[j,'N'])) # (access-use)/access
    tmp[j,'gap2']<-mean(emplogit2(tmp[j,'P'],tmp[j,'N'])-emplogit2(tmp[j,'Pu'],tmp[j,'N'])) # emplogit difference of access-use
    
    # not sure what this does
    tmp2<-HH[wh,]
    
    acc<-acc_mean<-agap<-ugap<-rep(NA,sum(wh))
    P<-d[wh]
    Pu<-e[wh]
    N<-a[wh]
    for(k in 1:nrow(tmp2)){ # for each household in cluster? --abv
      y<-tmp2$year
      h<-tmp2$n.individuals.that.slept.in.surveyed.hhs
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

# renaming --abv
HH1<-output
data<-HH1


# Cleanup: remove flawed points, print summary messages, save ------------------------------------------------------------

# what is complete.cases? --abv 
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
