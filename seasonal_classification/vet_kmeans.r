library(data.table)
library(plotly)


set.seed(12)
root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")

continent <- "africa"
this_cov <- "tsi_rainfall_vector_abundance_rescaled"
k <- 8

pal <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#000000")

main_dir <- file.path(base_dir, continent)

load(file.path(main_dir, "kmeans", paste0("k_out_", this_cov, "_", k, ".rdata")))
rotations <- fread(file.path(main_dir, paste0("svd_rotations_", this_cov, ".csv")))
traces <- fread(file.path(main_dir, "kmeans", paste0("random_trace_", this_cov, "_", k, ".csv")))

traces <- merge(traces, rotations, by="id", all.x=T)

centers <- data.table(k_out$centers)
centers[, cluster:="centroid"]
for_plot <- rbind(unique(traces[, list(cluster, X1, X2, X3)]), centers)
for_plot[, cluster:=as.factor(cluster)]

# ggplot(traces, aes(x=X1, y=X2, color=cluster)) +
#   geom_point(size=3, alpha=0.75) + 
#   scale_color_manual(values=pal) +
#   theme(legend.position = "none")


plot_ly(for_plot, x = ~X1, y = ~X2, marker=list(color="black"), z = ~X3, opacity=0.8) %>%
  add_markers()

plot_ly(for_plot, x = ~X1, y = ~X2, color = ~cluster, colors=pal, z = ~X3, opacity=0.8) %>%
  add_markers()
