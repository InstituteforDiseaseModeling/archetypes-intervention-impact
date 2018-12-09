library(data.table)
library(plotly)


set.seed(12)
root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")

continent <- "africa"
this_cov <- "tsi_rainfall_vector_abundance_rescaled"
k <- 8

main_dir <- file.path(base_dir, continent)

load(file.path(main_dir, "kmeans", paste0("k_out_", this_cov, "_", k, ".rdata")))
rotations <- fread(file.path(main_dir, paste0("svd_rotations_", this_cov, ".csv")))
traces <- fread(file.path(main_dir, "kmeans", paste0("random_trace_", this_cov, "_", k, ".csv")))

traces <- merge(traces, rotations, by="id", all.x=T)

plot_ly(traces, x = ~X1, y = ~X2, z = ~X3, color = ~cluster) %>%
  add_markers()


