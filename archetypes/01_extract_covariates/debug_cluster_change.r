

library(gdistance)
library(data.table)
library(stringr)
library(stats)
library(Hmisc)
library(ggplot2)
library(rasterVis)

rm(list=ls())

theme_set(theme_minimal(base_size = 16))
root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/")

overwrite <- T
rescale <- F
out_subdir <- "original_megatrends"
this_continent <- "africa"

out_dir <- file.path(base_dir, "results", out_subdir)
guide <- fread(file.path(out_dir, "instructions.csv"))

cov_details <- fread("map_covariates.csv")

cov_dir <- file.path(base_dir, "covariates", unique(guide$cov_directory), this_continent)
these_covs <- unlist(strsplit(guide[continent==this_continent]$covariates, "/"))

## Rasters themselves -----------------------------------------------------------------------------------------------------------------------

# this_cov_raster <- raster(file.path(cov_dir, "rainfall_chirps/rainfall_chirps_month_01.tif"))
# older_raster <- raster("~/Downloads/rainfall_chirps_month_01.tif")




## Extraction -----------------------------------------------------------------------------------------------------------------------


compare_vals <- lapply(these_covs, function(this_cov){
  
  print(this_cov)
  this_cov_dir <- file.path(cov_dir, this_cov)
  this_cov_details <- cov_details[cov==this_cov]
  
  this_cov_vals <- lapply(unlist(strsplit(this_cov_details$variables, "/")), function(this_variable){
    
    print(this_variable)
    
    this_cov_raster <- raster(file.path(this_cov_dir, paste0(this_cov, "_", this_cov_details$variable_label, "_", this_variable, ".tif")))
    
    version_one <- as.matrix(this_cov_raster)
    version_one_dt <- data.table(cov = this_cov,
                                 variable_name = this_cov_details$variable_label,
                                 variable_val = this_variable,
                                 version_one_id = which(!is.na(version_one)),
                                 version_one_value = version_one[!is.na(version_one)])
    
    version_two <- getValues(this_cov_raster)
    version_two_dt <- data.table(cov = this_cov,
                                 variable_name = this_cov_details$variable_label,
                                 variable_val = this_variable,
                                 version_two_id = which(!is.na(version_two)),
                                 version_two_value = version_two[!is.na(version_two)])
    
    for_id_map <- rep(NA, ncell(this_cov_raster))
    for_id_map[version_one_dt$version_one_id] <- version_one_dt$version_one_id
    for_id_map <- matrix(for_id_map, nrow=nrow(this_cov_raster))
    for_id_map <- base::t(for_id_map)
    id_map <- data.table(version_two_id = which(!is.na(for_id_map)),
                         version_one_id = for_id_map[!is.na(for_id_map)])

    full_comparison <- merge(version_one_dt, id_map, by="version_one_id", all=T)
    full_comparison <- merge(full_comparison, version_two_dt, by=c("cov", "variable_name", "variable_val", "version_two_id"), all=T)
    full_comparison[, diff:= abs(version_one_value - version_two_value)]
    
    return(full_comparison)
  })
  this_cov_vals <- rbindlist(this_cov_vals)
  return(this_cov_vals)

})

# keep only those pixels with values for all covariates
non_null_ids <- lapply(compare_vals, function(df){
  return(unique(df$version_one_id))
})
shared_ids <- Reduce(intersect, non_null_ids)
compare_vals <- rbindlist(compare_vals)
compare_vals <- compare_vals[version_one_id %in% shared_ids]

## SVD -----------------------------------------------------------------------------------------------------------------------

compare_vals_version_one <- compare_vals[, list(cov, variable_name, variable_val, id=version_one_id, value=version_one_value)]
print("reshaping and filling nulls")
svd_wide_datatable_version_one <- dcast(compare_vals_version_one, cov + variable_name + variable_val ~ id)
print("running svd")
svd_out_version_one <- svd(svd_wide_datatable_version_one[, 4:ncol(svd_wide_datatable_version_one)])

compare_vals_version_two <- compare_vals[, list(cov, variable_name, variable_val, id=version_two_id, value=version_two_value)]
print("reshaping and filling nulls")
svd_wide_datatable_version_two <- dcast(compare_vals_version_two, cov + variable_name + variable_val ~ id)
print("running svd")
svd_out_version_two <- svd(svd_wide_datatable_version_two[, 4:ncol(svd_wide_datatable_version_two)])

abs(svd_out_version_one$u[,1:3] - svd_out_version_two$u[, 1:3])
# same

all_vals_version_one <- melt(svd_wide_datatable_version_one, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
all_vals_version_one[, id:=as.integer(as.character(id))]

all_vals_version_two <- melt(svd_wide_datatable_version_two, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
all_vals_version_two[, id:=as.integer(as.character(id))]

## K-means -----------------------------------------------------------------------------------------------------------------------
nvecs <- 3
sing_vecs_version_one <- svd_out_version_one$u[, 1:nvecs]
sing_vecs_version_two <- svd_out_version_two$u[, 1:nvecs]

## multiply by original matrix to get rotations
print("rotating")
rotation_version_one <- data.frame(t(t(sing_vecs_version_one)%*%as.matrix(svd_wide_datatable_version_one[,4:ncol(svd_wide_datatable_version_one)])))
rotation_version_one$id <- as.integer(rownames(rotation_version_one))
rotation_version_one <- data.table(rotation_version_one)

print("rotating")
rotation_version_two <- data.frame(t(t(sing_vecs_version_two)%*%as.matrix(svd_wide_datatable_version_two[,4:ncol(svd_wide_datatable_version_two)])))
rotation_version_two$id <- as.integer(rownames(rotation_version_two))
rotation_version_two <- data.table(rotation_version_two)

id_map <- unique(compare_vals[, list(version_one_id, version_two_id)])

rotation_differences <- merge(rotation_version_one[, list(version_one_id=id, 
                                                          X1_v1=X1,
                                                          X2_v1=X2,
                                                          X3_v1=X3)],
                              id_map, by="version_one_id", all=T)

rotation_differences <- merge(rotation_differences,
                              rotation_version_two[, list(version_two_id=id, 
                                                          X1_v2=X1,
                                                          X2_v2=X2,
                                                          X3_v2=X3)],
                              by="version_two_id", all=T)

rotation_differences[, X1_diff:=abs(X1_v1-X1_v2)]
rotation_differences[, X2_diff:=abs(X2_v1-X2_v2)]
rotation_differences[, X3_diff:=abs(X3_v1-X3_v2)]

nclust <- 6
temp_raster <- raster(file.path(cov_dir, "mask.tif"))
these_colors <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A")


print(paste("finding clusters for k of", nclust))
load("original_seed.rdata")
.Random.seed <- seed_to_use
k_out_version_one <- kmeans(rotation_version_one[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
rotation_version_one[, cluster:= k_out_version_one$cluster]

# k_out_version_one[[1]] <- kmeans(rotation_version_one[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
# seed_to_use <- .Random.seed
# k_out_version_one[[2]] <- kmeans(rotation_version_one[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
# .Random.seed <- seed_to_use
# k_out_version_one[[3]] <- kmeans(rotation_version_one[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
# .Random.seed <- seed_to_use
# k_out_version_one[[4]] <- kmeans(rotation_version_one[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
# 
# k_outs_version_one <- lapply(k_out_version_one, function(this_k_out){
#   
#   this_rotation <- copy(rotation_version_one)
#   this_rotation[, cluster:=this_k_out$cluster]
#   
#   this_raster <- rep(NA, ncell(temp_raster))
#   this_raster[this_rotation$id] <- this_rotation$cluster
#   this_raster <- matrix(this_raster, nrow=nrow(temp_raster))
#   this_raster <- raster(this_raster, template=temp_raster)
#   this_raster <- ratify(this_raster)
#   
#   return(this_raster)
# })
# levelplot(stack(k_outs_version_one), att="ID", col.regions=these_colors,
#           xlab=NULL, ylab=NULL, scales=list(draw=F),
#           main = "", colorkey=F, margin=F)


set.seed(206)
k_out_version_two <- kmeans(rotation_version_two[, 1:nvecs], centers=nclust, algorithm = "MacQueen", iter.max=100)
rotation_version_two[, cluster:= k_out_version_two$cluster]


k_means_differences <- merge(rotation_version_one[, list(version_one_id=id, 
                                                         version_one_cluster=cluster)],
                             id_map, by="version_one_id", all=T)
k_means_differences <- merge(k_means_differences,
                             rotation_version_two[, list(version_two_id=id, 
                                                         version_two_cluster=cluster)],
                             by="version_two_id", all=T)

# here's the issue: 
k_means_differences[version_one_cluster != version_two_cluster]

# prove it with plots: 
print("creating new raster")
cluster_raster_version_one <- rep(NA, ncell(temp_raster))
cluster_raster_version_one[rotation_version_one$id] <- rotation_version_one$cluster
cluster_raster_version_one <- matrix(cluster_raster_version_one, nrow=nrow(temp_raster))
cluster_raster_version_one <- raster(cluster_raster_version_one, template=temp_raster)
cluster_raster_version_one <- ratify(cluster_raster_version_one)

cluster_raster_version_two <- rep(NA, ncell(temp_raster))
cluster_raster_version_two[rotation_version_two$id] <- rotation_version_two$cluster
cluster_raster_version_two <- matrix(cluster_raster_version_two, nrow=nrow(temp_raster), byrow=T)
cluster_raster_version_two <- raster(cluster_raster_version_two, template=temp_raster)
cluster_raster_version_two <- ratify(cluster_raster_version_two)

levelplot(stack(cluster_raster_version_one, cluster_raster_version_two), att="ID", col.regions=these_colors,
                      xlab=NULL, ylab=NULL, scales=list(draw=F),
                      main = "", colorkey=F, margin=F)


levelplot(cluster_raster_version_two, att="ID", col.regions=these_colors,
          xlab=NULL, ylab=NULL, scales=list(draw=F),
          main = "", colorkey=F, margin=F)


#whatt? compare to raster output from last time: 
cluster_raster <- raster(file.path(out_dir, this_continent, "02_kmeans", paste0("map_", nclust, "_cluster", ".tif")))
cluster_raster <- ratify(cluster_raster)
levelplot(cluster_raster, att="ID", col.regions=these_colors,
                      xlab=NULL, ylab=NULL, scales=list(draw=F),
                      main = "", colorkey=F, margin=F) 

# yup, still different. trace backwards
load(file.path(out_dir, this_continent, "02_kmeans", paste0("k_out_", nclust, "_cluster", ".rdata")))
k_out$centers

# try rotation
orig_rotations <- fread(file.path(out_dir, this_continent, "02_kmeans", "svd_rotations.csv"))


# svd 
load(file.path(out_dir, this_continent, "01_svd", "svd_output.rdata"))
sing_vecs_orig <- svd_out$u[, 1:nvecs]

test_new_svd <- svd(svd_wide_datatable[, 4:ncol(svd_wide_datatable)])

# try reordering data

all_vals_orig <- melt(svd_wide_datatable, id.vars=c("cov", "variable_name", "variable_val"), variable.name="id", value.name="cov_val")
all_vals_orig[, id:=as.integer(as.character(id))]
all_vals_orig[, variable_val:= str_pad(variable_val, 2, side="left", pad="0")]
all_vals_orig <- all_vals_orig[order(cov, variable_name, variable_val)]

svd_wide_orig_reordered <- dcast(all_vals_orig, cov + variable_name + variable_val ~ id)
test_reordered_svd <- svd(svd_wide_orig_reordered[, 4:ncol(svd_wide_orig_reordered)])

test <- all_vals_version_two[, list(cov, variable_name, variable_val, id, cov_val_test=cov_val)]
big_test <- merge(test, all_vals_orig, by=c("cov", "variable_name", "variable_val", "id"), all=T)
big_test[, diff:=abs(cov_val_test- cov_val)]


