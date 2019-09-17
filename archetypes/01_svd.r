## -----------------------------------------------------------------------------------------------------------------
# Seasonality Classification
# 01_svd.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# June 2018
# 
# Extending the work described in the link below, this script takes the tsi and rainfall covariates 
# extracted in extract_raster_values and runs svd on the joint dataset. 
# 
# https://paper.dropbox.com/doc/Cluster-MAP-pixels-by-seasonality-zga4UM1DnBx8pc11rStOS
## -----------------------------------------------------------------------------------------------------------------------

library(gdistance)
library(data.table)
library(stringr)
library(stats)
library(Hmisc)
library(ggplot2)

rm(list=ls())

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification/")
continents <- c("africa", "asia", "americas")
cov_details <- fread("clustering_covariates.csv")
cov_details <- cov_details[variable!="year"]
overwrite <- F
rescale <- T
theme_set(theme_minimal(base_size = 16))

for (continent in continents){
  print(paste("running svd on", continent))
  main_dir <- file.path(base_dir, continent)
  
  these_covs <- cov_details[continents %like% continent]
  full_label <- paste(these_covs$cov, collapse="_")
  if (rescale==T & "rainfall" %in% these_covs$cov){
    full_label <- paste0(full_label, "_rescaled")
  }
  
  # load datasets, merge
  all_vals_fname <- file.path(main_dir, paste0(full_label, "_vals.csv"))
  if (file.exists(all_vals_fname) & overwrite==F){
    print("loading extracted values")
    all_vals <- fread(all_vals_fname)
  }else{
    print("appending datasets")
    all_vals <- lapply(these_covs$cov, function(cov_name){
      vals <- fread(file.path(main_dir, paste0(cov_name, "_vals.csv")))
      if (rescale==T & cov_name=="rainfall"){
        
        # cap outliers in rainfall data; rescale to [0,1]
        vals[value>400, value:=400]
        vals[, value:=value/max(value)]
      }
      return(vals)
    })
    
    # keep only those pixels with values for all covariates
    non_null_ids <- lapply(all_vals, function(df){
      return(unique(df$id))
    })
    shared_ids <- Reduce(intersect, non_null_ids)
    all_vals <- rbindlist(all_vals)
    all_vals <- all_vals[id %in% shared_ids]
    write.csv(all_vals, file=all_vals_fname, row.names=F)
  }
  
  # plot distribution of values
  png(file=file.path(main_dir, paste0(full_label, "_distributions.png")))
  
  distplot <- ggplot(all_vals, aes(x=value)) + 
              geom_density(aes(color=cov, fill=cov), alpha=0.5) +
              facet_grid(cov~.) +
              theme_minimal() + 
              theme(legend.position = "none") + 
              labs(title="Distribution of SVD variables",
                   x="Rescaled Covariate Value",
                   y="Density") 
  print(distplot)
  
  graphics.off()
  
  # svd
  svd_out_fname <- file.path(main_dir, paste0(full_label, "_svd.rdata"))
  if (file.exists(svd_out_fname) & overwrite==F){
    print("loading svd outputs")
    load(svd_out_fname)
  }else{
    print("reshaping and filling nulls")
    for_svd <- dcast(all_vals, cov + variable_name + variable_val ~ id)
    print("running svd")
    svd_out <- svd(for_svd[, 4:ncol(for_svd)])
    save(svd_out, file=svd_out_fname)
  }
  
  print("plotting")
  init_variance <- svd_out$d^2/sum(svd_out$d^2)
  variance <- data.table(continent=capitalize(continent),
                         vector=1:length(init_variance), 
                         variance_explained=init_variance)
  
  pdf(file=file.path(main_dir, paste0(full_label, "_svd.pdf")))
  varplot <- ggplot(variance[vector<=5], aes(x=vector, y=variance_explained)) +
    geom_line(size=2) +
    geom_point(size=5) +
    theme(legend.position = "none") +
    labs(x="Singular Vector", 
         y="Variance Explained"
         #title=paste("Variance Explained by Singular Vectors", capitalize(continent))
         )
  print(varplot)
  
  graphics.off()
  print(varplot)
}



