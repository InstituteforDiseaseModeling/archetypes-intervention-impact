
# STEP 5: PREDICT

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-64 --logging gs://map_data_z/users/amelia/logs --input-recursive in_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_rerun_sam joint_dir=gs://map_data_z/users/amelia/itn_cube/input_data_archive/  cov_dir=gs://map_data_z/cubes_5km func_dir=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud --input CODE=gs://map_data_z/users/amelia/itn_cube/code/run_on_gcloud/05_itn_prediction.r  --output-recursive out_dir=gs://map_data_z/users/amelia/itn_cube/results/20190606_rerun_sam/05_predictions --command 'Rscript ${CODE}'

# from access deviation and use gap, predict coverage and map all surfaces

package_load <- function(package_list){
  # package installation/loading
  new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  lapply(package_list, library, character.only=T)
}

package_load(c("zoo","raster", "data.table", "rgdal", "INLA", "RColorBrewer",
               "VGAM", "maptools", "rgeos", "png", "sp", "doParallel", "rasterVis", "ggplot2"))

if(Sys.getenv("in_dir")=="") {
  in_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/"
  joint_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data"
  cov_dir <- "/Volumes/GoogleDrive/Team Drives/cubes/5km incomplete/"
  out_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190507_sam_withseeds/05_predictions"
  func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
} else {
  in_dir <- Sys.getenv("in_dir") 
  joint_dir <- Sys.getenv("joint_dir")
  cov_dir <- Sys.getenv("cov_dir") # cubes_5km
  out_dir <- Sys.getenv("out_dir")
  func_dir <- Sys.getenv("func_dir") # code directory for function scripts
}

# load shapefiles and templates for plotting
World<-readShapePoly(file.path(joint_dir, 'for_plots_and_tables/MAP_templates_Africa', 'afri_word.shp'))
World <- gSimplify(World, tol=0.1, topologyPreserve=TRUE)
Africa<-readShapePoly(file.path(joint_dir, 'for_plots_and_tables/MAP_templates_Africa', 'Africa.shp'))
Africa <- gSimplify(Africa, tol=0.1, topologyPreserve=TRUE)
Water<-readShapePoly(file.path(joint_dir, 'for_plots_and_tables/MAP_templates_Africa', 'africa_water.shp'))
Water <- gSimplify(Water, tol=0.1, topologyPreserve=TRUE)
img <- readPNG(file.path(joint_dir, 'for_plots_and_tables/map-logo_bg.png'))
print("loaded shapefiles and templates for plotting")

# load access and gap estimates
list.files(in_dir)
load(file.path(in_dir, "03_access_deviation.Rdata"))
mod.pred.acc=mod.pred
theta.acc<-theta
mesh.acc<-mesh

load(file.path(in_dir, "04_use_gap.Rdata"))
mod.pred.use=mod.pred
theta.use<-theta
mesh.use<-mesh

#overwrite like-named functions in INLAfunctions.r
source(file.path(func_dir, "05_itn_prediction_functions.r"))


# HUGE Prediction Loop ----------------------------------------------------------------------------####################### 
time_points<-2000:2016
static.covs<-get.pred.covariates.static()

for(xxx in 1:length(time_points)){
  
  dics<-cors<-c()
  pred_year=time_points[xxx]
  print(paste("predicting for", pred_year))
  pred_month1="01"
  pred_month2="02"
  pred_month3="03"
  pred_month4="04"
  pred_month5="05"
  pred_month6="06"
  pred_month7="07"
  pred_month8="08"
  pred_month9="09"
  pred_month10="10"
  pred_month11="11"
  pred_month12="12"
  pred_year_cov=pred_year
  if(time_points[xxx]<=2000) pred_year_cov=2001 # covariates are not full in 2000
  if(time_points[xxx]>=2013) pred_year_cov=2013 # covariates are not full in 2000
  
  pred_mat_y<-get.pred.covariates.yearonly(pred_year_cov)
  
  ####################################################################################################################
  #extract covariates
  print("extracting covariates")
  pred_mat1_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month1)
  pred_mat2_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month2)
  pred_mat3_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month3)
  pred_mat4_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month4)
  pred_mat5_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month5)
  pred_mat6_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month6)
  pred_mat7_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month7)
  pred_mat8_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month8)
  pred_mat9_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month9)
  pred_mat10_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month10)
  pred_mat11_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month11)
  pred_mat12_s<-get.pred.covariates.dynamic(pred_year_cov,pred_month12)
  
  
  pred_mat1<-cbind(static.covs,pred_mat_y,pred_mat1_s)
  pred_mat2<-cbind(static.covs,pred_mat_y,pred_mat2_s)
  pred_mat3<-cbind(static.covs,pred_mat_y,pred_mat3_s)
  pred_mat4<-cbind(static.covs,pred_mat_y,pred_mat4_s)
  pred_mat5<-cbind(static.covs,pred_mat_y,pred_mat5_s)
  pred_mat6<-cbind(static.covs,pred_mat_y,pred_mat6_s)
  pred_mat7<-cbind(static.covs,pred_mat_y,pred_mat7_s)
  pred_mat8<-cbind(static.covs,pred_mat_y,pred_mat8_s)
  pred_mat9<-cbind(static.covs,pred_mat_y,pred_mat9_s)
  pred_mat10<-cbind(static.covs,pred_mat_y,pred_mat10_s)
  pred_mat11<-cbind(static.covs,pred_mat_y,pred_mat11_s)
  pred_mat12<-cbind(static.covs,pred_mat_y,pred_mat12_s)
  print("covariates extracted")
  
  ####################################################################################################################
  #standardise and transform covariates
  #for(i in 1:12){	
  #	for(j in 1:length(trans)) {
  #		eval(parse(text=paste0("pred_mat",i,"[,",j,"]<-t[[trans[",j,"]]](pred_mat",i,"[,",j,"])")))
  #		eval(parse(text=paste0("pred_mat",i,"[,",j,"]<-(pred_mat",i,"[,",j,"]-ms[",j,"])/sds[",j,"]")))
  #	}
  #}
  ####################################################################################################################
  print("loading prediction locations and inla results")
  #get pred locations
  pred_locs<-get.pred.locs()
  ####################################################################################################################
  #get field
  #########################################
  pred_year=time_points[xxx]	     ########
  pred_year_field=pred_year
  if(pred_year>=2015) pred_year_field=2015 # adjust prediction year for field
  
  #########################################
  INLA:::inla.dynload.workaround() 
  
  pred_year1=pred_year
  Apreds<-inla.spde.make.A(mesh.acc, loc=cbind(pred_locs[,'x'],pred_locs[,'y'],pred_locs[,'z']),group=rep(pred_year_field,nrow(pred_locs)),group.mesh=mesh1d)
  
  A.pred_dynamic1=Apreds
  
  
  print("prediction locations and inla results loaded")
  ####################################################################################################################
  #### get access mean from stock and flow
  print("loading stock and flow means")
  cn<-raster(file.path(joint_dir, 'general/african_cn5km_2013_no_disputes.tif')) #load raster
  NAvalue(cn)=-9999
  pred_val<-getValues(cn)#get values again
  w<-is.na(pred_val) #find NAs again
  index<-1:length(w) 
  index<-index[!w]
  gauls<-cn[index]
  
  
  load(file.path(joint_dir, 'stock_and_flow/prop0prop1.rData'))
  POPULATIONS<-read.csv(file.path(joint_dir, 'general/National_Config_Data.csv'))
  Country.list<-Cout
  Country.gaul<-Cout
  for(i in 1:length(Country.gaul)){
    Country.gaul[i]<-POPULATIONS[as.character(POPULATIONS[,'ISO3'])==as.character(Country.list[i]),'GAUL_Code']
  }	
  times_ind<-seq(0.0,18,0.25)+2000
  Country.gaul = as.numeric(as.vector(Country.gaul))
  accdate1<-as.numeric(as.yearmon(paste(time_points[xxx],1,sep='-')))
  accdate2<-as.numeric(as.yearmon(paste(time_points[xxx],2,sep='-')))
  accdate3<-as.numeric(as.yearmon(paste(time_points[xxx],3,sep='-')))
  accdate4<-as.numeric(as.yearmon(paste(time_points[xxx],4,sep='-')))
  accdate5<-as.numeric(as.yearmon(paste(time_points[xxx],5,sep='-')))
  accdate6<-as.numeric(as.yearmon(paste(time_points[xxx],6,sep='-')))
  accdate7<-as.numeric(as.yearmon(paste(time_points[xxx],7,sep='-')))
  accdate8<-as.numeric(as.yearmon(paste(time_points[xxx],8,sep='-')))
  accdate9<-as.numeric(as.yearmon(paste(time_points[xxx],9,sep='-')))
  accdate10<-as.numeric(as.yearmon(paste(time_points[xxx],10,sep='-')))
  accdate11<-as.numeric(as.yearmon(paste(time_points[xxx],11,sep='-')))
  accdate12<-as.numeric(as.yearmon(paste(time_points[xxx],12,sep='-')))
  
  
  acc1_p<-rep(NA,length(gauls))
  acc2_p<-rep(NA,length(gauls))
  acc3_p<-rep(NA,length(gauls))
  acc4_p<-rep(NA,length(gauls))
  acc5_p<-rep(NA,length(gauls))
  acc6_p<-rep(NA,length(gauls))
  acc7_p<-rep(NA,length(gauls))
  acc8_p<-rep(NA,length(gauls))
  acc9_p<-rep(NA,length(gauls))
  acc10_p<-rep(NA,length(gauls))
  acc11_p<-rep(NA,length(gauls))
  acc12_p<-rep(NA,length(gauls))	
  
  print("stock and flow means loaded")
  
  print("calculating country-level access")
  for(i in 1:length(Country.list)){
    print(Country.list[i])
    p<-out[[i]]
    p0<-p[,,1] # get p0 parameters
    p1<-p[,,2] # get p1 parameters	
    func1<-func0<-list()
    # make these into continuous functions
    for(f in 1:ncol(p0)){
      eval(parse(text=paste0('func0[[',f,']]<-approxfun(times_ind,p0[,',f,'])')))
    }
    for(f in 1:ncol(p0)){
      eval(parse(text=paste0('func1[[',f,']]<-approxfun(times_ind,p1[,',f,'])')))
    }
    # load household size distribution
    hh<-read.csv(file.path(joint_dir, 'stock_and_flow/HHsize.csv'))
    KEY=read.csv(file.path(joint_dir, 'stock_and_flow/hh_surveys_key.csv'))
    cn_nm<-as.character(KEY[KEY$Name%in%Country.list[i],'Svy.Name'])
    hh_val<-hh[hh$HHSurvey%in%cn_nm,]
    if(nrow(hh_val)!=0){
      print("USING COUNTRY SUBSET FOR HOUSEHOLD SIZE DISTRIBUTION")
      hh_val<-colSums(hh_val[,3:ncol(hh_val)])/sum(colSums(hh_val[,3:ncol(hh_val)]))
    } else {
      print("USING ENTIRE HSIZE FOR HOUSEHOLD SIZE DISTRIBUTION")
      hh_val<-colSums(hh[,3:ncol(hh)])/sum(colSums(hh[,3:ncol(hh)]))
    }
    hh_size<-hh_val[1:10]
    hh_size[10]<-sum(hh_val[10:length(hh_val)])
    hh<-hh_size
    names(hh)<-1:10
    
    
    #### calculate access
    
    acc1<-calc_access(accdate1,func0,func1,hh)
    acc2<-calc_access(accdate2,func0,func1,hh)
    acc3<-calc_access(accdate3,func0,func1,hh)
    acc4<-calc_access(accdate4,func0,func1,hh)
    acc5<-calc_access(accdate5,func0,func1,hh)
    acc6<-calc_access(accdate6,func0,func1,hh)
    acc7<-calc_access(accdate7,func0,func1,hh)
    acc8<-calc_access(accdate8,func0,func1,hh)
    acc9<-calc_access(accdate9,func0,func1,hh)
    acc10<-calc_access(accdate10,func0,func1,hh)
    acc11<-calc_access(accdate11,func0,func1,hh)
    acc12<-calc_access(accdate12,func0,func1,hh)
    
    cc<-gauls==Country.gaul[i]
    acc1_p[cc]=emplogit(acc1,1000)	
    acc2_p[cc]=emplogit(acc2,1000)
    acc3_p[cc]=emplogit(acc3,1000)
    acc4_p[cc]=emplogit(acc4,1000)
    acc5_p[cc]=emplogit(acc5,1000)
    acc6_p[cc]=emplogit(acc6,1000)
    acc7_p[cc]=emplogit(acc7,1000)
    acc8_p[cc]=emplogit(acc8,1000)
    acc9_p[cc]=emplogit(acc9,1000)
    acc10_p[cc]=emplogit(acc10,1000)
    acc11_p[cc]=emplogit(acc11,1000)
    acc12_p[cc]=emplogit(acc12,1000)
    
  }
  print("country-level access calculated")
  ##### put it all together
  
  print("predicting access for new locations")
  Beta<-c(mod.pred.acc$summary.fixed[,1])
  Beta_cov<-Beta[2:length(Beta)]
  
  
  
  fixed_effects1<-Beta[1] + as.numeric(Beta_cov%*%t(pred_mat1))	
  fixed_effects2<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat2))	
  fixed_effects3<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat3))	
  fixed_effects4<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat4))	
  fixed_effects5<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat5))	
  fixed_effects6<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat6))	
  fixed_effects7<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat7))	
  fixed_effects8<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat8))	
  fixed_effects9<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat9))	
  fixed_effects10<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat10))	
  fixed_effects11<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat11))	
  fixed_effects12<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat12))	
  
  
  lp1=fixed_effects1 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp2=fixed_effects2 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp3=fixed_effects3 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp4=fixed_effects4 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp5=fixed_effects5 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp6=fixed_effects6 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp7=fixed_effects7 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp8=fixed_effects8 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp9=fixed_effects9 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp10=fixed_effects10 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp11=fixed_effects11 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  lp12=fixed_effects12 + drop(A.pred_dynamic1%*%mod.pred.acc$summary.random$field$mean) 
  
  
  lp1<-Inv.IHS(lp1,theta.acc)
  lp2<-Inv.IHS(lp2,theta.acc)
  lp3<-Inv.IHS(lp3,theta.acc)
  lp4<-Inv.IHS(lp4,theta.acc)
  lp5<-Inv.IHS(lp5,theta.acc)
  lp6<-Inv.IHS(lp6,theta.acc)
  lp7<-Inv.IHS(lp7,theta.acc)
  lp8<-Inv.IHS(lp8,theta.acc)
  lp9<-Inv.IHS(lp9,theta.acc)
  lp10<-Inv.IHS(lp10,theta.acc)
  lp11<-Inv.IHS(lp11,theta.acc)
  lp12<-Inv.IHS(lp12,theta.acc)
  
  
  lp1_lp = acc1_p+lp1
  lp2_lp = acc2_p+lp2
  lp3_lp = acc3_p+lp3
  lp4_lp = acc4_p+lp4
  lp5_lp = acc5_p+lp5
  lp6_lp = acc6_p+lp6
  lp7_lp = acc7_p+lp7
  lp8_lp = acc8_p+lp8
  lp9_lp = acc9_p+lp9
  lp10_lp = acc10_p+lp10
  lp11_lp = acc11_p+lp11
  lp12_lp = acc12_p+lp12
  
  acc_mean<-(acc1_p + acc2_p + acc3_p +acc4_p+acc5_p + acc6_p + acc7_p + acc8_p + acc9_p + acc10_p + acc11_p +acc12_p)/12
  
  lp = cbind(lp1,lp2,lp3,lp4,lp5,lp6,lp7,lp8,lp9,lp10,lp11,lp12) ## access deviation
  
  lp_lp = cbind(lp1_lp,lp2_lp,lp3_lp,lp4_lp,lp5_lp,lp6_lp,lp7_lp,lp8_lp,lp9_lp,lp10_lp,lp11_lp,lp12_lp) # access
  
  
  lp<-(lp)
  lp<-rowMeans(lp)
  
  lp_lp<-plogis(lp_lp)
  lp_lp<-rowMeans(lp_lp)
  
  
  colfunc <- colorRampPalette(c("blue","cyan","yellow","orange","red"))
  colfunc2 <- colorRampPalette(c("blue","cyan","yellow","orange","red"),bias=2)
  
  P<-cn
  P[!w]<-(lp_lp)
  P[!is.na(cn) & is.na(P)]=0
  
  P2<-cn
  P2[!w]<-(lp)
  P2[!is.na(cn) & is.na(P2)]=0
  
  P3<-cn
  P3[!w]<-plogis(acc_mean)
  P3[!is.na(cn) & is.na(P3)]=0
  
  print("access, access deviation, and mean access calculated, saving rasters")
  
  writeRaster(P, file.path(out_dir, paste0('ITN_',time_points[xxx],'.ACC.tif')),NAflag=-9999,overwrite=TRUE)
  writeRaster(P2, file.path(out_dir, paste0('ITN_',time_points[xxx],'.DEV.tif')),NAflag=-9999,overwrite=TRUE)
  writeRaster(P3, file.path(out_dir, paste0('ITN_',time_points[xxx],'.MEAN.tif')),NAflag=-9999,overwrite=TRUE)
  
  #############################################################################################################################################
  #now predict use
  
  print("predicting use")
  pred_year1=pred_year
  Apreds<-inla.spde.make.A(mesh.use, loc=cbind(pred_locs[,'x'],pred_locs[,'y'],pred_locs[,'z']),group=rep(pred_year,nrow(pred_locs)),group.mesh=mesh1d)
  
  A.pred_dynamic1=Apreds
  
  
  Beta<-c(mod.pred.use$summary.fixed[,1])
  Beta_cov<-Beta[2:length(Beta)]
  
  
  fixed_effects1<-Beta[1] + as.numeric(Beta_cov%*%t(pred_mat1))	
  fixed_effects2<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat2))	
  fixed_effects3<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat3))	
  fixed_effects4<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat4))	
  fixed_effects5<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat5))	
  fixed_effects6<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat6))	
  fixed_effects7<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat7))	
  fixed_effects8<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat8))	
  fixed_effects9<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat9))	
  fixed_effects10<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat10))	
  fixed_effects11<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat11))	
  fixed_effects12<-Beta[1] +as.numeric(Beta_cov%*%t(pred_mat12))	
  
  
  gap1=fixed_effects1 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap2=fixed_effects2 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap3=fixed_effects3 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap4=fixed_effects4 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap5=fixed_effects5 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap6=fixed_effects6 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap7=fixed_effects7 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap8=fixed_effects8 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap9=fixed_effects9 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap10=fixed_effects10 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap11=fixed_effects11 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  gap12=fixed_effects12 + drop(A.pred_dynamic1%*%mod.pred.use$summary.random$field$mean) 
  
  
  gap1<-Inv.IHS(gap1,theta.use)
  gap2<-Inv.IHS(gap2,theta.use)
  gap3<-Inv.IHS(gap3,theta.use)
  gap4<-Inv.IHS(gap4,theta.use)
  gap5<-Inv.IHS(gap5,theta.use)
  gap6<-Inv.IHS(gap6,theta.use)
  gap7<-Inv.IHS(gap7,theta.use)
  gap8<-Inv.IHS(gap8,theta.use)
  gap9<-Inv.IHS(gap9,theta.use)
  gap10<-Inv.IHS(gap10,theta.use)
  gap11<-Inv.IHS(gap11,theta.use)
  gap12<-Inv.IHS(gap12,theta.use)
  
  use1<-(lp1_lp-gap1)
  use2<-(lp2_lp-gap2)
  use3<-(lp3_lp-gap3)
  use4<-(lp4_lp-gap4)
  use5<-(lp5_lp-gap5)
  use6<-(lp6_lp-gap6)
  use7<-(lp7_lp-gap7)
  use8<-(lp8_lp-gap8)
  use9<-(lp9_lp-gap9)
  use10<-(lp10_lp-gap10)
  use11<-(lp11_lp-gap11)
  use12<-(lp12_lp-gap12)
  
  lp_gap = cbind(gap1,gap2,gap3,gap4,gap5,gap6,gap7,gap8,gap9,gap10,gap11,gap12) # access
  lp_gap<-rowMeans(lp_gap)
  
  
  lp_use = cbind(use1,use2,use3,use4,use5,use6,use7,use8,use9,use10,use11,use12) # access
  lp_use<-plogis(lp_use)
  lp_use<-rowMeans(lp_use)
  # predict
  P4<-cn
  P4[!w]<-(lp_gap)
  P4[!is.na(cn) & is.na(P4)]=0
  
  P5<-cn
  P5[!w]<-(lp_use)
  P5[!is.na(cn) & is.na(P5)]=0
  
  print("use predicted, saving rasters")
  
  writeRaster(P4, file.path(out_dir, paste0('ITN_',time_points[xxx],'.GAP.tif')),NAflag=-9999,overwrite=TRUE)
  writeRaster(P5, file.path(out_dir, paste0('ITN_',time_points[xxx],'.USE.tif')),NAflag=-9999,overwrite=TRUE)
  
  
}


## todo: new script here


# Plotting ----------------------------------------------------------------------------####################### 

print("plotting results")

time_points<-2000:2016
colfunc <- colorRampPalette(c("red","orange","yellow",'green'))

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

gainr<-c(-2,1) # relative gain?
data<-read.csv(file.path(in_dir, '01_database.csv'))

# Plot access ----
filename<-file.path(out_dir,'access.pdf')
file.remove(filename)
pdf(filename,width=18.61, height=10.34,onefile=T)

for(i in time_points){
  P=raster(file.path(out_dir, paste0('ITN_',i,'.ACC.tif')),NAflag=-9999)
  #	P2=raster(file.path(out_dir, paste0('ITN_',i,'.DEV.tif')),NAflag=-9999)
  P3=raster(file.path(out_dir, paste0('ITN_',i,'.MEAN.tif')),NAflag=-9999)
  #	P4=raster(file.path(out_dir, paste0('ITN_',i,'.GAP.tif')),NAflag=-9999)
  P5=raster(file.path(out_dir, paste0('ITN_',i,'.USE.tif')),NAflag=-9999)
  relg<-relative_gain(P5,P)
  relg_decision<-relg
  relg_decision[relg_decision>0]=1
  relg_decision[relg_decision<0]=-1
  relg_decision[relg_decision==0]=0
  
  P[1:5]=1;P[6:10]=0
  P3[1:5]=1;P3[6:10]=0
  P5[1:5]=1;P5[6:10]=0
  #	P2[1:5]=devr[1];P2[6:10]=devr[2]
  #	P4[1:5]=gapr[1];P4[6:10]=gapr[2]
  relg[1:5]=gainr[1];relg[6:10]=gainr[2]
  relg[relg<gainr[1]]=gainr[1]
  P6<-(P-P5)/P
  P6=1-P6
  P6[P6<1e-2]=0
  
  par(mfrow=c(2,4),mar=c(1,1,1,1))
  Plot_raster(P,1,"ITN access")
  #	Plot_raster(P2,1,"ITN access deviation")
  Plot_raster(P3,1,"ITN access mean")
  #	Plot_raster(P4,1,"ITN use gap")
  Plot_raster(P5,1,"ITN use")
  Plot_raster(P6,1,"Off Acces")
  Plot_raster(relg,0.3,"ITN relative gain")
  Plot_raster(relg_decision,1,"ITN decision gain")
  Plot_points(data,i,'use','Raw use data')
  Plot_points(data,i,'access','Raw access data')
  
  
}
graphics.off()

# Plot relative gain (?) ----
#use
use=P3
access=P	
r=use-access
r[r<=0]=0
l=r+P

#access
use=P3
access=P	
r=use-access
l2=1+r
l2[l2>=1]=1

#use
use=P3 #(proportion population sleeping under a net)
access=P #(proportion population with access to an ITN in their household)	 
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

par(mfrow=c(1,2))

plot(relativegain,col=jet.colors(1000),main='relative gain')

relativegain[relativegain>0]=1
relativegain[relativegain<0]=-1
relativegain[relativegain==0]=0

plot(relativegain,main='relative gain decision')


# Generate country-level estimates ----
ncores <- detectCores()
print(paste("--> Machine has", ncores, "cores available"))
registerDoParallel(ncores-2)

tables<-foreach(i=time_points) %dopar% {	
  print("loading just-created access and use rasters")
  P=raster(file.path(out_dir, paste0('ITN_',i,'.ACC.tif')),NAflag=-9999)
  P2=raster(file.path(out_dir, paste0('ITN_',i,'.MEAN.tif')),NAflag=-9999)
  P3=raster(file.path(out_dir, paste0('ITN_',i,'.USE.tif')),NAflag=-9999)
  print("just-created rasters successfully loaded")
  print("loading population raster")
  print(file.path(cov_dir, paste0('AfriPop/',min(2015,i),'.total.population.tif')))
  # ABV EDIT: there is not 2015 raster, so we have to cap the pop loading at 2015.
  pop<-raster(file.path(cov_dir, paste0('AfriPop/', min(2015,i),'.total.population.tif')),NAflag=-9999)
  print("population raster loaded")
  P4=(P-P2)/P2 # relative access deviation
  P5=(P3-P)/P # relative use gap
  P6=(P-P2) # absolute access deviation
  P7=(P3-P) # absolute use gap
  
  t1<-get.cn.estimates(P,pop)
  t2<-get.cn.estimates(P2,pop)
  t3<-get.cn.estimates(P3,pop)
  t4<-get.cn.estimates(P4,pop)
  t5<-get.cn.estimates(P5,pop)
  t6<-get.cn.estimates(P6,pop)
  t7<-get.cn.estimates(P7,pop)
  
  return(cbind(t1,t2,t3,t4,t5,t6,t7))
  
}	
print("end dopar")

# Plot use as a function of access ----

print("plotting summary.pdf")
filename<-file.path(out_dir, 'Summary.pdf')
file.remove(filename)
pdf(filename,width=11.59, height=10.06,onefile=T) 
i=2013
for(i in time_points){
  
  P=raster(file.path(out_dir, paste0('ITN_',i,'.ACC.tif')),NAflag=-9999)
  P2=raster(file.path(out_dir, paste0('ITN_',i,'.MEAN.tif')),NAflag=-9999)
  P5=raster(file.path(out_dir, paste0('ITN_',i,'.USE.tif')),NAflag=-9999)
  
  P[1:5]=1;P[6:10]=0
  P2[1:5]=1;P2[6:10]=0
  P5[1:5]=1;P5[6:10]=0
  
  P6<-(P-P5)/P
  P6=1-P6
  P6[P6<1e-2]=0
  
  par(mfrow=c(2,2),mar=c(1,1,1,1))
  Plot_raster(P,1,"ITN access")
  Plot_raster(P2,1,"ITN Stock and flow mean")
  Plot_raster(P5,1,"ITN use")
  Plot_raster(P6,1,"Use as function of Access")
  
}
graphics.off()
print("summary.pdf plotted")

# Final adjustments and plotting (ask Sam)----

print("running final adjustments and plotting")
# this routine adjust use numbers to zero in bad years and label files the same as before for backward compatibility
print("loading indicators")
indicators<-read.csv(file.path(joint_dir, 'stock_and_flow/indicators_access_qtr_new.csv'))

ind<-matrix(nrow=nrow(indicators),ncol=17)
for(i in 1:nrow(indicators)){
  ind[i,]<-tapply(as.numeric(indicators[i,2:ncol(indicators)]),c(rep(1:(72/4), each = 4),18),mean)[1:17]
}
rownames(ind)<-indicators[,1]
colnames(ind)<-seq(2000,2016,by=1)
indicators<-ind

# todo: does this POPULATIONS pull a different file than the one earlier in the script?
print("loading populations")
POPULATIONS<-read.csv(file.path(joint_dir, 'general/country_table_populations.csv')) # load table to match gaul codes to country names
names<-as.character(rownames(ind))
#for(i in 1:nrow(indicators)){
#	names[i]=as.character(POPULATIONS[as.character(POPULATIONS$NAME)==names[i],'COUNTRY_ID']) # get 3 letter country codes
#}
print("generating indicator matrix")
rownames(ind)<-as.character(names)
threshold=0.02
indmat<-matrix(data=0,nrow=nrow(indicators),ncol=ncol(indicators))
for(i in 1:nrow(indicators)){
  for(j in 1:ncol(indicators)){
    if(indicators[i,j]<threshold){
      indmat[i,j]=1
    }
    if(indicators[i,j]>threshold){
      break
    }
  }
}
colnames(indmat)<-seq(2000,2016,by=1)
rownames(indmat)<-names

times<-seq(2000,2016,by=1)

fnames_for_stack <- paste(out_dir, '/ITN_',times,'.USE.tif',sep="")
print("contents of stack directory:")
print(out_dir)
list.files(out_dir)
print("loading raster stack")
print(fnames_for_stack)

st<-stack(paste(out_dir, '/ITN_',times,'.USE.tif',sep=""))
print("plotting original use stack")

pal <- brewer.pal(8, "RdYlGn")
breaks <- seq(0, 1, length.out = length(pal)+1)

pdf(file.path(out_dir, "samdir_zdir_compare.pdf"), width=7, height=10)
p <- levelplot(st,
               par.settings=rasterTheme(region=pal), at=breaks,
               xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F,
               main="Initial Use Rasters"
)
print(p)


stv<-getValues(st)
cn<-raster(file.path(joint_dir, 'general/african_cn5km_2013_no_disputes.tif'))
cnv<-getValues(cn)
for(i in 1:nrow(indmat)){
  gaul<-POPULATIONS[as.character(POPULATIONS$NAME)==rownames(ind)[i],'GAUL_CODE']
  tmp<-cn
  wh<-colnames(indmat)[indmat[i,]==1]
  wh<-as.numeric(wh)-1999
  stv[cnv==gaul,wh]<-0
}

for(i in 1:nlayers(st)){
  st[[i]]<-setValues(st[[i]],stv[,i])
}

p <- levelplot(st,
               par.settings=rasterTheme(region=pal), at=breaks,
               xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F,
               main="Use Rasters for Z Drive"
)
print(p)
graphics.off()

print("saving for z dir rasters")

for(i in 1:nlayers(st)){
  writeRaster(st[[i]],file.path(out_dir, paste("for_z_", times[i],'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
  
}

# create zero years
st0<-st[[1]]
NAvalue(st0)<--9999
st0[!is.na(st0)]=0
writeRaster(st0, file.path(out_dir, paste("for_z_", 1999,'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
writeRaster(st0, file.path(out_dir, paste("for_z_", 1998,'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
writeRaster(st0, file.path(out_dir, paste("for_z_", 1997,'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
writeRaster(st0, file.path(out_dir, paste("for_z_", 1996,'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
writeRaster(st0, file.path(out_dir, paste("for_z_", 1995,'.ITN.use.yearavg.new.adj.tif',sep="")),NAflag=-9999,overwrite=TRUE)
print("for z dir rasters saved")

