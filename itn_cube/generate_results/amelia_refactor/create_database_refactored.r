
## Amelia's first pass at parsing "Create_database_ITN access....r"

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
# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-16 --logging gs://map_data_z/users/amelia/logs --input-recursive input_dir=gs://map_data_z/users/amelia/itn_cube/create_database/input joint_dir=gs://map_data_z/users/amelia/itn_cube/joint_data --input func_fname=gs://map_data_z/users/amelia/itn_cube/code/create_database_functions.r CODE=gs://map_data_z/users/amelia/itn_cube/code/create_database_refactored.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/create_database/output --command 'Rscript ${CODE}'

# Data loading, household-level access/use stats  ------------------------------------------------------------

if(Sys.getenv("input_dir")=="") {
  # gcsfuse drives are just too laggy
  # joint_dir <- "/Users/bertozzivill/Desktop/zdrive_mount/users/amelia/itn_cube/joint_data"
  # input_dir <- "/Users/bertozzivill/Desktop/zdrive_mount/users/amelia/itn_cube/create_database/input"
  # output_dir <- "/Users/bertozzivill/Desktop/zdrive_mount/users/amelia/itn_cube/create_database/output"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/input"
  output_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/output"
  func_fname <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/create_database_functions.r"
} else {
  joint_dir <- Sys.getenv("joint_dir")
  input_dir <- Sys.getenv("input_dir")
  output_dir <- Sys.getenv("output_dir")
  func_fname <- Sys.getenv("func_fname")
}

source(func_fname)
out_fname <- file.path(output_dir, 'ITN_final_clean_access_18March2019_COMPARE.csv')

# p0 & p1-- from stock & flow, nat'l time series of p0=p(hh has 0 nets) and p1=avg # of nets
# 40 countries (list length), houshold size 1-10 (columns)
load(file.path(input_dir, "prop0prop1.rData")) # contains "out" (list of s$f time series) and "Cout" (country names)
stock_and_flow_isos<-Cout # Cout is a vector of ISOs in prop0prop1.rData
names(out) <- stock_and_flow_isos


# load household data and survey-to-country key, keep only those in country list
survey_to_country_key <- read.csv(file.path(input_dir, "KEY_080817.csv"))  #TODO: is this used?
HH<-fread(file.path(input_dir, "ALL_HH_Data_20112017.csv")) # todo: come back and delete cols we don't need. also rename this
HH<-HH[ISO3 %in% stock_and_flow_isos]
# TODO: remove missing values up here instead of at the end

# Remove households with size zero or NA
HH <- HH[!is.na(n.individuals.that.slept.in.surveyed.hhs) & n.individuals.that.slept.in.surveyed.hhs>0]

# Find unique surveys
unique_surveys <- unique(HH$Survey.hh)

# find access (# with a net available) and use (# sleeping under net) per household 
times<-seq(2000, 2018,0.25) # quarter-year intervals

#  update HH with columns for the above values
# todo: check bounding of d against methods in document
# a: already present as n.individuals.that.slept.in.surveyed.hhs
# b: already present as n.ITN.per.hh
# d: turn into n.covered.by.nets
HH[, n.covered.by.nets:=pmin(n.ITN.per.hh*2, n.individuals.that.slept.in.surveyed.hhs)]
# e: already present as n.individuals.that.slept.under.ITN

# # test locally
# orig_unique_surveys <- copy(unique_surveys)
# unique_surveys <- c("SN2014DHS")


# Main loop: calculating access/gap for each household cluster  ------------------------------------------------------------

ncores <- detectCores()
print(paste("--> Machine has", ncores, "cores available"))
registerDoParallel(ncores-2)
output<-foreach(i=1:length(unique_surveys),.combine=rbind) %dopar% { #survey loop
  
  svy<-unique_surveys[i] # store survey
  
  this_survey_data=HH[Survey.hh==svy,] # keep only household data for the survey in question
  un<-unique(this_survey_data$Cluster.hh) # get unique cluster IDs for that household survey
  
  ## find distribution of household sizes in this survey as measured by the number of people sleeping under a net the night before
  this_survey_data[, sample.prop:=sample.w/sum(sample.w)]
  this_survey_data[, capped.n.sleeping.in.hhs:=pmin(n.individuals.that.slept.in.surveyed.hhs, 10)]
  household_props <- this_survey_data[, list(hh.size.prop=sum(sample.prop)), by=list(capped.n.sleeping.in.hhs)]
  household_props <- household_props[order(capped.n.sleeping.in.hhs)] # used to be called "hh"
  
  svy_country <- unique(this_survey_data[Survey.hh==svy]$ISO3)
  s_and_f_probs<-out[[svy_country]]
  # both p0 and p1 are 69 x 10 datasets with row as time point and column as household size
  p0<-s_and_f_probs[,,1] # get p0: probability of no nets in household
  p1<-s_and_f_probs[,,2] # get p1: average # of nets in household

  # Sam defines a set of functions that linearly interpolate between these quarterly values, 
  # but then goes on to only evaluate it at rounded years. Can we use the month information 
  # from the dataset to do a better job?
  
  # make these into continuous functions via linear interpolation
  func0 <- lapply(1:ncol(p0), function(col){
    approxfun(times, p0[,col])})
  func1 <- lapply(1:ncol(p1), function(col){
    approxfun(times, p1[,col])})
  
  # define access for this survey-year based on household props
  # todo: refactor calc_access_matrix function
  access_dt <- lapply(unique(this_survey_data$year), function(this_year){
    this_survey_access <- calc_access_matrix(this_year, func0, func1, household_props$hh.size.prop)
    this_access_dt <- household_props[, list(year=this_year, capped.n.sleeping.in.hhs, 
                                        from_matrix_access=this_survey_access[1:10],
                                        access_matrix_mean=this_survey_access[11])]
    
  })
  access_dt <-rbindlist(access_dt)
  
  this_survey_data <- merge(this_survey_data, access_dt, by=c("year", "capped.n.sleeping.in.hhs"), all.x=T)
  
  
  # aggregate to cluster level
  summary_by_cluster <- this_survey_data[, list(Survey=svy,
                                                lat=mean(latitude),
                                                lon=mean(longitude),
                                                P=sum(n.covered.by.nets),
                                                N=sum(n.individuals.that.slept.in.surveyed.hhs),
                                                year=mean(year),
                                                Pu=sum(n.individuals.that.slept.under.ITN),
                                                T=sum(n.ITN.per.hh),
                                                gap3=mean(1-n.individuals.that.slept.under.ITN/n.covered.by.nets, na.rm=T),
                                                Amean=mean(from_matrix_access),
                                                Tmean=mean(access_matrix_mean)
                                                ),
                                         by=list(Cluster.hh)]
  summary_by_cluster[, gap:=( (P/N)-(Pu/N) ) / (P/N)] # (access-use)/access
  summary_by_cluster[, gap2:=emplogit2(P,N) - emplogit2(Pu, N)] # emplogit difference of access-use
  summary_by_cluster <- summary_by_cluster[order(Cluster.hh)]
  summary_by_cluster[, Cluster.hh:=NULL]
  
  return(summary_by_cluster)
}

# renaming
data<-data.frame(output)


# Cleanup: remove flawed points, print summary messages, save ------------------------------------------------------------

#print(paste('**OUTPUT MESSAGE** remove points with no cooridnates: there are - ',nrow(data[!complete.cases(data),])))
data<-data[complete.cases(data),] # NOTE: do you want to be getting rid of cases where the only NA's are in the "gap" values?

print(paste('**OUTPUT MESSAGE** remove points with 0 lat 0 lon: there are - ',nrow(data[data$lat==0 & data$lon==0,])))
data<-data[data$lat!=0 & data$lon!=0,]

# Check for invalid points, and attempt to reposition them
cn<-raster(file.path(joint_dir, 'african_cn5km_2013_no_disputes.tif')) # master country layer
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
