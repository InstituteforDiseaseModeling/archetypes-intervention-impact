rm(list=ls())

library(data.table)
library(ggplot2)
theme_set(theme_minimal(base_size=14))

# main_dir <- "C:/Users/abertozzivilla/repos/map-larval-habitat/covariate_prep/"
out_dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/Mozambique/incidence_calibration/transmission_population_comparison/"
main_dir <- out_dir


data <- fread(paste0(main_dir, "prev_data_2_10.csv"))
map <- fread(paste0(main_dir, "map_prev.csv"))

setnames(map, "map_prev_2015", "prev")
map <- map[grid_cell %in% unique(data$grid_cell)]
data <- data[, list(type="Data", grid_cell, mid_x, mid_y, catchment, prev, pop)]
data <- merge(data, unique(map[, list(grid_cell, map_pop)]), by="grid_cell", all=T)

all <- rbind(data, map)
all <- all[, list(grid_cell, type, mid_x, mid_y, catchment, prev, 
                  census_pop=pop, log_census_pop=log10(pop), afripop=map_pop, log_afripop=log10(map_pop))]
means <- all[, list(prev=mean(prev), census_pop=mean(census_pop), afripop=mean(afripop),
                    log_census_pop=log10(mean(census_pop)), log_afripop=log10(mean(afripop))),
             by=list(type, catchment)]

for (val in c("census_pop", "afripop", "log_census_pop", "log_afripop")){
  
  label <- ifelse(val=="census_pop", "Census Pop", "Afripop")
  if (val %like% "log"){
    label <- paste0("log10(", label, ")")
  }
  print(label)
  
  pdf(paste0(out_dir, "prev_", val, ".pdf"), height=9, width=12)
  plot <- ggplot(all, aes_string(x="prev", y=val, color="catchment")) +
          geom_point(alpha=0.25) +
          geom_point(data=means,aes_string(x="prev", y=val), size=3) +
          facet_wrap(~type) +
          scale_color_brewer(type="qual", palette="Paired") +
          labs(title=paste("Magude Prevalence by Catchment, Prevalence, and", label),
               x="Prevalence Ages 2-10",
               y=label)
  if (val %like% "log"){
    plot <- plot + ylim(0, max(all$log_census_pop))
  }else{
    plot <- plot + ylim(0, max(all$census_pop))
  }
  
  print(plot)
  graphics.off()
  
  
}






