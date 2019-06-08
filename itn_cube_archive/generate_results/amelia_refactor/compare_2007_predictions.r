
library(data.table)

rm(list=ls())

load("/Volumes/GoogleDrive/My Drive/itn_cube/results/20190507_sam_withseeds/05_acc_image_2007.Rdata")
new_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/results/20190521_replicate_prediction/"
national_raster_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/joint_data/african_cn5km_2013_no_disputes.tif"
func_dir <- "/Users/bertozzivill/repos/malaria-atlas-project/itn_cube/generate_results/amelia_refactor/"
source(file.path(func_dir, "check_file_similarity.r"))


new_predict_dir <- file.path(new_dir, "05_predictions")
# access, access deviation, use gap

# national access (not pixel)
new_nat_access <- fread(file.path(new_predict_dir, "national_access_2007_for_comparison.csv"))

test <- new_nat_access[iso3=="ZWE"]
test <- cbind(test, c(acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10, acc11, acc12))
test[, diff:=nat_access-V2]

# national access (pixel)
nat_access <- list(acc1_p, acc2_p, acc3_p, acc4_p, acc5_p, acc6_p, acc7_p, acc8_p, acc9_p, acc10_p, acc11_p, acc12_p)
acc_dev_predictions <- fread(file.path(new_predict_dir, "/acc_dev_2007_for_comparison.csv"))
acc_dev_predictions <- acc_dev_predictions[order(cellnumber)]

for (this_month in 1:12){
  print(paste("comparing nat access for month", this_month))
  new_subset <- acc_dev_predictions[month==this_month]
  old_subset <- data.table(emplogit_nat_access=nat_access[[this_month]])
  check_sameness(old_subset, new_subset)
}

# access deviation by month

fixed_effects <- list(fixed_effects1, fixed_effects2, fixed_effects3, fixed_effects4, fixed_effects5, fixed_effects6,
                      fixed_effects7, fixed_effects8, fixed_effects9, fixed_effects10, fixed_effects11, fixed_effects12)

full_dev <- list(lp1, lp2, lp3, lp4, lp5, lp6, lp7, lp8, lp9, lp10, lp11, lp12)

IHS <- function(x, theta){  # function to IHS transform
  return(asinh(theta * x)/theta) 
}

for (this_month in 1:12){
  print(paste("comparing access stats for month", this_month))
  new_subset <- acc_dev_predictions[month==this_month]
  old_subset <- data.table(fixed=fixed_effects[[this_month]],
                           access_deviation=full_dev[[this_month]])
  old_subset[, random:=IHS(access_deviation, theta.acc)-fixed]
  check_sameness(old_subset, new_subset, sameness_cutoff=7e-05)
}

# summary values for raster

old_nat <- cbind(acc1_p, acc2_p, acc3_p, acc4_p, acc5_p, acc6_p, acc7_p, acc8_p, acc9_p, acc10_p, acc11_p, acc12_p)
old_dev <- cbind(lp1,lp2,lp3,lp4,lp5,lp6,lp7,lp8,lp9,lp10,lp11,lp12) ## access deviation
old_acc <- cbind(lp1_lp,lp2_lp,lp3_lp,lp4_lp,lp5_lp,lp6_lp,lp7_lp,lp8_lp,lp9_lp,lp10_lp,lp11_lp,lp12_lp) # access

transform_old <- data.table(nat_access=rowMeans(plogis(old_nat)),
                            access_deviation=rowMeans(plogis(old_dev)),
                            access=rowMeans(plogis(old_acc)))


transform_acc_dev_predictions <- acc_dev_predictions[, list(iso3, year, month, cellnumber, 
                                                  emplogit_access=emplogit_access,
                                                  access=plogis(emplogit_access),
                                                  access_deviation=plogis(access_deviation),
                                                  nat_access=plogis(emplogit_nat_access))]

summary_access <- transform_acc_dev_predictions[, list(nat_access=mean(nat_access, na.rm=F),
                                             access_deviation=mean(access_deviation, na.rm=F),
                                             access=mean(access, na.rm=F)),
                                      by=list(iso3, year, cellnumber)]
summary_access <- summary_access[order(cellnumber)]

check_sameness(transform_old, summary_access, sameness_cutoff = 7e-05)

# using weirdly transformed variables

old_sam_odd <- data.table(nat_access=acc_mean,
                          access_deviation=lp,
                          access=lp_lp)

sam_odd <- acc_dev_predictions[, list(nat_access=mean(emplogit_nat_access),
                                      access_deviation=mean(access_deviation)), by=list(iso3, year, cellnumber)]
new_sam_odd <- sam_odd[, list(iso3, year, cellnumber, nat_access,
                          access_deviation, 
                          access=summary_access$access)]

check_sameness(old_sam_odd, new_sam_odd, sameness_cutoff = 7e-05)

# tifs
national_raster <- raster(national_raster_dir)
NAvalue(national_raster) <- -9999

access_map <- copy(national_raster)
access_map[summary_access$cellnumber] <- summary_access$access
access_map[!is.na(national_raster) & is.na(access_map)] <- 0

access_diff <- access_map-P

deviation_map <- copy(national_raster)
deviation_map[new_sam_odd$cellnumber] <- new_sam_odd$access_deviation
deviation_map[!is.na(national_raster) & is.na(deviation_map)] <- 0
deviation_diff <- deviation_map-P2

nat_mean_map <- copy(national_raster)
nat_mean_map[new_sam_odd$cellnumber] <- plogis(new_sam_odd$nat_access)
nat_mean_map[!is.na(national_raster) & is.na(nat_mean_map)] <- 0
nat_mean_diff <- nat_mean_map-P3




# covariates

new_dynamic <- fread(file.path(new_dir, "02_dynamic_covariates/dynamic_2007.csv"))
old_dynamic <- list(pred_mat1_s, pred_mat2_s, pred_mat3_s, pred_mat4_s, pred_mat5_s,
                    pred_mat6_s, pred_mat7_s, pred_mat8_s, pred_mat9_s, pred_mat10_s,
                    pred_mat11_s, pred_mat12_s)

for (this_month in 1:12){
  print(paste("comparing for month", this_month))
  new_subset <- new_dynamic[month==this_month]
  old_subset <- data.table(old_dynamic[[this_month]])
  check_sameness(old_subset, new_subset)
}


pred_mat_y <- data.table(pred_mat_y)
new_annual <- fread(file.path(new_dir, "02_annual_covariates.csv"))
new_annual <- new_annual[year==2007]
check_sameness(pred_mat_y, new_annual)


static.covs <- data.table(static.covs)
new_static <- fread(file.path(new_dir, "02_static_covariates.csv"))
check_sameness(static.covs, new_static)