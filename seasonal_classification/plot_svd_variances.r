# quick(ish) and dirty: nicer SVD plots

library(data.table)
library(ggplot2)
library(Hmisc)

main_dir <- file.path(Sys.getenv("HOME"), "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/seasonal_classification")

continents <- c("asia", "americas", "africa")

# variances <- data.table()

variances <- lapply(continents, function(continent){
  print(continent)
  fname <- ifelse(continent=="africa", "tsi_rainfall_vector_abundance_svd.rdata", "tsi_rainfall_svd.rdata")
  load(file.path(main_dir, continent, fname))
  init_variance <- svd_out$d^2/sum(svd_out$d^2)
  variance <- data.table(continent=capitalize(continent),
                         vector=1:length(init_variance), 
                         variance_explained=init_variance)
  return(variance)
})

variances <- rbindlist(variances)

pdf("/Users/bertozzivill/Desktop/svd_curves.pdf", height=3, width=7)
ggplot(variances[vector<=10], aes(x=vector, y=variance_explained, color=continent)) +
  geom_line(size=1) +
  geom_point(shape=1, size=2) +
  scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
  theme_minimal() +
  facet_grid(~continent) +
  theme(legend.position = "none") +
  labs(x="Singular Vector", 
       y="Variance Explained",
       title="SVD Outputs by Continent: Variance Explained by Singular Vectors")

graphics.off()
