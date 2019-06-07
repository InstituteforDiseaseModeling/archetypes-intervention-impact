
## STEP ONE: CREATE DATABASE

## sam: takes Harry's data, does aggregating, find national level means at quarterly times 

rm(list=ls())
tic <- Sys.time()

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster","VGAM", "doParallel", "data.table"))

# current dsub:
# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-16 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/input_data_archive --input func_fname=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud/01_create_database_functions.r CODE=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud/01_create_database.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_rerun_sam --command 'Rscript ${CODE}'

# Data loading, household-level access/use stats  ------------------------------------------------------------

if(Sys.getenv("input_dir")=="") {
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/input"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/output"
  func_fname <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/create_database_functions.r"
} else {
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir")
  func_fname <- Sys.getenv("func_fname")
}

source(func_fname)

out_fname <- file.path(output_dir, '01_database.csv')

# p0 & p1-- from stock & flow, nat'l time series of p0=p(hh has >0 nets) and p1=avg # of nets
# 40 countries (list length), houshold size 1-10 (columns)
load(file.path(input_dir, "stock_and_flow/prop0prop1.rData"))
KEY <- read.csv(file.path(input_dir, "stock_and_flow/props_key.csv"))

# Cout is a vector of ISOs in prop0prop1.rData
Country.list<-Cout

# load household data, keep only those in country list
HH<-read.csv(file.path(input_dir, "database/ALL_HH_Data_20112017.csv"),stringsAsFactors=F)
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
# todo: give these better names, put them in a dataset? check bounding of d against methods in document
times<-seq(0.0,18,0.25)+2000 # 2000 to 2018, in quarter-year intervals
a=HH$n.individuals.that.slept.in.surveyed.hhs # number in household
b=HH$n.ITN.per.hh # number of nets
d=b*2 # number covered by nets
d[d>a]=a[d>a] # bounding to not exceed the number in the house
e=HH$n.individuals.that.slept.under.ITN # use count

# # test locally
# orig_Surveys <- copy(Surveys)
# Surveys <- Surveys[2:3]

# Main loop: calculating access/gap for each household cluster  ------------------------------------------------------------

ncores <- detectCores()
print(paste("--> Machine has", ncores, "cores available"))
registerDoParallel(ncores-2)
output<-foreach(i=1:length(Surveys),.combine=rbind) %dopar% { #survey loop
  svy<-Surveys[i] # store survey
  
  cc<-HH[HH$Survey.hh==svy,'ISO3'] # get country
  p<-out[[which(Country.list==cc[1])]]
  # both p0 and p1 are 69 x 10 datasets with row as time point and column as household size
  p0<-p[,,1] # get p0: probability of at least 1 net in household
  p1<-p[,,2] # get p1: average # of nets in household
  
  func0<-list()
  func1<-list()
  # make these into continuous functions (via linear interpolation--abv)
  for(f in 1:ncol(p0)){
    eval(parse(text=paste0('func0[[',f,']]<-approxfun(times,p0[,',f,'])')))
  }
  for(f in 1:ncol(p0)){ # should this be p1? may not matter if they're the same length
    eval(parse(text=paste0('func1[[',f,']]<-approxfun(times,p1[,',f,'])')))
  }
  
  tmp_hh=HH[HH$Survey.hh==svy,] # keep only household data for the survey in question
  un<-unique(tmp_hh$Cluster.hh) # get unique cluster IDs for that household survey
  
  ### house hold size as measured by the number of people sleeping under a net the night before
  ## (find distribution of household sizes in this survey -- abv)
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
  # todo: transform this into the direct creation of a data.frame from list of cluster ids
  # p,n,t: data, means: means from stock and flow
  fields<-c("Survey" ,"lat" , "lon",'year','P','N','Pu','T','Amean','Tmean','gap','gap2','gap3')
  tmp<-matrix(NA,nrow=length(un),ncol=length(fields))
  colnames(tmp)<-fields
  tmp<-as.data.frame(tmp)
  
  for(j in 1:length(un)){ #for each cluster
    # binary indicator, same length as HH, of "is this the survey cluster you want? (todo: just subset HH?)
    wh<-HH$Survey.hh==svy & HH$Cluster.hh==un[j]
    
    tmp[j,'Survey']<-svy
    tmp[j,'lat']<-mean(HH[wh,'latitude'],na.rm=TRUE)
    tmp[j,'lon']<-mean(HH[wh,'longitude'],na.rm=TRUE)
    
    tmp[j,'P']<-sum(d[wh]) # number with access to nets in this cluster
    tmp[j,'N']<-sum(a[wh]) # number in households in cluster
    tmp[j,'year']<-mean(HH$year[wh])
    tmp[j,'Pu']<-sum(e[wh]) # number using nets in this cluster
    tmp[j,'T']<-sum(b[wh]) # number of nets in this cluster
    tmp[j,'gap']<-mean(((tmp[j,'P']/tmp[j,'N'])-(tmp[j,'Pu']/tmp[j,'N']))/(tmp[j,'P']/tmp[j,'N'])) # (access-use)/access
    tmp[j,'gap2']<-mean(emplogit2(tmp[j,'P'],tmp[j,'N'])-emplogit2(tmp[j,'Pu'],tmp[j,'N'])) # emplogit difference of access-use
    
    # NOW you subset HH??
    tmp2<-HH[wh,]
    
    # initialize a bunch of empty vectors with as many entries as there is household data in HH for that cluster
    acc<-acc_mean<-agap<-ugap<-rep(NA,sum(wh))
    P<-d[wh] # number with access to nets
    Pu<-e[wh] # number using nets
    N<-a[wh] # people in households
    for(k in 1:nrow(tmp2)){ # for each household in cluster, why are you doing this to me
      y<-tmp2$year
      h<-tmp2$n.individuals.that.slept.in.surveyed.hhs
      h[h>10]=10
      access<-calc_access_matrix(y[k],func0,func1,hh) # unclear what exactly calc_access_matrix does
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

# renaming
HH1<-output
data<-HH1


# Cleanup: remove flawed points, print summary messages, save ------------------------------------------------------------

#print(paste('**OUTPUT MESSAGE** remove points with no cooridnates: there are - ',nrow(data[!complete.cases(data),])))
data<-data[complete.cases(data),] # NOTE: do you want to be getting rid of cases where the only NA's are in the "gap" values?

print(paste('**OUTPUT MESSAGE** remove points with 0 lat 0 lon: there are - ',nrow(data[data$lat==0 & data$lon==0,])))
data<-data[data$lat!=0 & data$lon!=0,]

# Check for invalid points, and attempt to reposition them
cn<-raster(file.path(input_dir, 'general/african_cn5km_2013_no_disputes.tif')) # master country layer
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

print(paste("--> Writing to", out_fname))
write.csv(data, out_fname, row.names=FALSE)

toc <- Sys.time()
elapsed <- toc-tic

print(paste("--> Time Elapsed: ", elapsed, units(elapsed)))