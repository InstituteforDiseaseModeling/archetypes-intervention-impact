library(data.table)
library(ggplot2)

## compare outcomes of different climate files on PfPR 2-10 annually


rm(list=ls())
theme_set(theme_minimal(base_size = 18))

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/compare_climate")


fname <- "MAP_ERA5_Climate_2000_2016_Inset_Burnin.csv"

all_data <- lapply(list.files(main_dir), function(fname){
  label <- sub("MAP_(.+)_Inset_Burnin\\.csv", "\\1", fname)
  data <- fread(file.path(main_dir, fname))
  data[,climate_type:=label]
  return(data)
})

all_data <- rbindlist(all_data)
all_data[, mean_prev:=mean(final_prev), by=list(climate_type, Site_Name, x_Temporary_Larval_Habitat, day)]
all_data[, year:=round(day/365)+1999]

mean_data <- unique(all_data[, list(climate_type, Site_Name, year, x_Temporary_Larval_Habitat, mean_prev)])
x_temps <- unique(mean_data$x_Temporary_Larval_Habitat)
mean_data[, x_temp:=round(x_Temporary_Larval_Habitat, 4)]
mean_data[, climate_type:=factor(climate_type, levels=c("IDM_Climate_2016", "ERA5_Climate_2016",
                                                        "ERA5_Climate_2000_2016", "ERA5_Climate_2000_2016_Landtemp"))]

for (site in unique(mean_data$Site_Name)){
  print(ggplot(mean_data[Site_Name==site], aes(x=year, y=mean_prev, color=climate_type)) + 
    geom_line() +
    facet_wrap(x_temp ~ .) +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    labs(title=site,
         x="Year",
         y="Mean PfPR 2-10")
    )
}

subset <- mean_data[x_temp>0.1 & x_temp<10]

ggplot(subset[x_temp==0.631], aes(x=year, y=mean_prev, color=climate_type)) + 
  geom_line(size=1) +
  facet_grid(. ~ Site_Name) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="Year",
       y="Mean PfPR 2-10")



final_results <- mean_data[year==2014]
ggplot(final_results, aes(x=x_Temporary_Larval_Habitat, y=mean_prev, color=climate_type)) +
  geom_point(size=1.5) +
  scale_x_continuous(trans='log10') + 
  facet_wrap(Site_Name ~ .) + 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="x_Temp",
       y="Mean PfPR 2-10")










