###############################################################################################################
## V2_00_calculate_intervention_classes.r
## Amelia Bertozzi-Villa
## April 2020
## 
## Using the most recent intervention surfaces, divide Africa into high/medium/low ITN/ACT categories
## and yes/no IRS categories
##############################################################################################################

library(data.table)
library(raster)
library(rasterVis)

rm(list=ls())

year_to_use <- 2019

analysis_subdir <- "20200426_int_history"
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")
out_dir <- file.path(base_dir, analysis_subdir, "input")


itn_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20200418_BMGF_ITN_C1.00_R1.00_V2/04_predictions"
act_dir <- "~/Desktop/covid_mitigation_project/effective_treatment_with_an_antimalarial3"
irs_dir <- "~/Desktop/covid_mitigation_project/IRS/"

itn_layer <- raster(file.path(itn_dir, paste0("ITN_", year_to_use, "_use.tif")))
act_layer <- raster(file.path(act_dir, paste0(year_to_use, ".effective.treatment.tif")))
irs_layer <- raster(file.path(irs_dir, paste0(year_to_use, ".IRS.tif")))

# mask to ITN layer, it's the most restricted
itn_layer[itn_layer==0] <- NA

irs_layer <- raster::mask(irs_layer, itn_layer)
act_layer <- raster::mask(act_layer, itn_layer)
act_layer[is.na(act_layer) & itn_layer>0] <- 0 #fill in stray NA's

# find quantiles
quant_cutoffs <- c(0.33, 0.66)
itn_quantiles <- quantile(itn_layer, probs=quant_cutoffs)
act_quantiles <- quantile(act_layer, probs=quant_cutoffs)

# reclassify

make_rcl <- function(these_quantiles){
  rcl <- matrix(data=c(0, these_quantiles[[1]], 1,
                       these_quantiles[[1]], these_quantiles[[2]], 2,
                       these_quantiles[[2]], 1, 3),
                nrow=3, byrow=T)
  return(rcl)
}

itn_categorical_layer <- reclassify(itn_layer, make_rcl(itn_quantiles), right=F)
act_categorical_layer <- reclassify(act_layer, make_rcl(act_quantiles), right=F)
irs_categorical_layer <- irs_layer>0


# make full 18-unit map
# initialize with raster of zeros
full_categorical <- copy(irs_categorical_layer)
full_categorical[full_categorical==1] <- 0

full_categorical[itn_categorical_layer==1 & act_categorical_layer==1 & irs_categorical_layer==0] <- 1
full_categorical[itn_categorical_layer==2 & act_categorical_layer==1 & irs_categorical_layer==0] <- 2
full_categorical[itn_categorical_layer==3 & act_categorical_layer==1 & irs_categorical_layer==0] <- 3

full_categorical[itn_categorical_layer==1 & act_categorical_layer==2 & irs_categorical_layer==0] <- 4
full_categorical[itn_categorical_layer==2 & act_categorical_layer==2 & irs_categorical_layer==0] <- 5
full_categorical[itn_categorical_layer==3 & act_categorical_layer==2 & irs_categorical_layer==0] <- 6

full_categorical[itn_categorical_layer==1 & act_categorical_layer==3 & irs_categorical_layer==0] <- 7
full_categorical[itn_categorical_layer==2 & act_categorical_layer==3 & irs_categorical_layer==0] <- 8
full_categorical[itn_categorical_layer==3 & act_categorical_layer==3 & irs_categorical_layer==0] <- 9

full_categorical[itn_categorical_layer==1 & act_categorical_layer==1 & irs_categorical_layer==1] <- 10
full_categorical[itn_categorical_layer==2 & act_categorical_layer==1 & irs_categorical_layer==1] <- 11
full_categorical[itn_categorical_layer==3 & act_categorical_layer==1 & irs_categorical_layer==1] <- 12

full_categorical[itn_categorical_layer==1 & act_categorical_layer==2 & irs_categorical_layer==1] <- 13
full_categorical[itn_categorical_layer==2 & act_categorical_layer==2 & irs_categorical_layer==1] <- 14
full_categorical[itn_categorical_layer==3 & act_categorical_layer==2 & irs_categorical_layer==1] <- 15

full_categorical[itn_categorical_layer==1 & act_categorical_layer==3 & irs_categorical_layer==1] <- 16
full_categorical[itn_categorical_layer==2 & act_categorical_layer==3 & irs_categorical_layer==1] <- 17
full_categorical[itn_categorical_layer==3 & act_categorical_layer==3 & irs_categorical_layer==1] <- 18

writeRaster(full_categorical, filename=file.path(out_dir, "intervention_history_categories.tif"), overwrite=T)

## find mean intervention coverage values for each category
category_labels <- data.table(int_category=1:18,
                              label = c(
                                "Low ITN, Low ACT, No IRS",
                                "Medium ITN, Low ACT, No IRS",
                                "High ITN, Low ACT, No IRS",
                                "Low ITN, Medium ACT, No IRS",
                                "Medium ITN, Medium ACT, No IRS",
                                "High ITN, Medium ACT, No IRS",
                                "Low ITN, High ACT, No IRS",
                                "Medium ITN, High ACT, No IRS",
                                "High ITN, High ACT, No IRS",
                                "Low ITN, Low ACT, Yes IRS",
                                "Medium ITN, Low ACT, Yes IRS",
                                "High ITN, Low ACT, Yes IRS",
                                "Low ITN, Medium ACT, Yes IRS",
                                "Medium ITN, Medium ACT, Yes IRS",
                                "High ITN, Medium ACT, Yes IRS",
                                "Low ITN, High ACT, Yes IRS",
                                "Medium ITN, High ACT, Yes IRS",
                                "High ITN, High ACT, Yes IRS"
                              ))

coverages_for_model <- rbindlist(lapply(1:18, function(this_category){
  print(this_category)
  
  data.table(int_category=this_category,
             itn=round(mean(itn_layer[full_categorical==this_category]),2),
             act=round(mean(act_layer[full_categorical==this_category]),2),
             irs= round(mean(irs_layer[full_categorical==this_category]),2),
             n_pixels=sum(full_categorical[full_categorical==this_category])
  )
}))

coverages_for_model <- merge(coverages_for_model, category_labels, by="int_category", all=T)

write.csv(coverages_for_model, file=file.path(out_dir, "intervention_history_categories.csv"), row.names = F)

for_pyscript <- melt(coverages_for_model, id.vars = "int_category", measure.vars = c("itn", "act", "irs"), 
                     variable.name="int", 
                     value.name="cov")
setnames(for_pyscript, "int_category", "int_id")
for_pyscript[int=="act", int:="al_cm"]
for_pyscript[, start_day:=0]
for_pyscript[, cov:=cov*100]
write.csv(for_pyscript, file=file.path(out_dir, "interventions.csv"), row.names = F)







