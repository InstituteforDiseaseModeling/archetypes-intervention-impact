library(data.table)
library(ggplot2)

rm(list=ls())

main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/prelim_itn_sweeps")

lookup_karen <- fread(file.path(main_dir, "lookup_table_karen_multi_int.csv"))
lookup_moine <- fread(file.path(main_dir, "lookup_table_moine_multi_int.csv"))

data_all <- fread(file.path(main_dir, "full_karen_data.csv"))

data_all <- data_all[year<max(year) & Run_Number==0 & year<50] 
data_all[, infected:=PfPR_by_Age_Bin*Average_Population_by_Age_Bin]


age_dist <- data_all[year %in% c(0, 10, 20, 30, 40, 50), list(pop=mean(Average_Population_by_Age_Bin)), by=list(age, year)]

# make equivalent age bins
even_ages <- age_dist[age<1, list(age=0, pop=sum(pop)), by="year"]

for (ageval in 1:109){
  print(ageval)
  if (ageval<15){
    new_data <- age_dist[age==ageval]
  }else if (ageval<30){
    new_data <- age_dist[age>=15 & age<=ageval]
    new_data[, pop:=pop/5]
  }else{
    new_data <- age_dist[age>=30 & age<=ageval]
    new_data[, pop:=pop/10]
  }
  
  even_ages <- rbind(even_ages, new_data)
  
}



age_dist[, perc:=pop/sum(pop), by="year"]
ggplot(even_ages, aes(x=age, y=pop)) +
  geom_point() + 
  facet_grid(~year) +
  labs(title="Age Distribution by Year",
       x="Age Bin",
       y="Population")
sum(age_dist[age<10 & age>=2]$pop)/sum(age_dist$pop)


ggplot(data_all, aes(x=year, y=Average_Population_by_Age_Bin)) +
  geom_line(aes(color=factor(x_Temporary_Larval_Habitat))) +
  facet_wrap(~age, scales="free_y") +
  theme(legend.position = "none")



# png(file.path(main_dir, "anthro_sweep.png"), width=3, height=2, units="in", res=100)
# ggplot(lookup_karen, aes(x=initial, y=final)) +
#   geom_line(size=1.5) +
#   # facet_grid(.~type) +
#   theme_minimal() +
#   labs(title="SE Asia ITN+ACT Impact",
#        x="Initial Prevalence",
#        y="Final Prevalence")
# graphics.off()

lookup_karen[, ITN_Coverage:=factor(ITN_Coverage)]
lookup_karen[, ACT_Coverage:= factor(ACT_Coverage)]

ggplot(lookup_karen[IRS_Coverage==0], aes(x=initial, y=final)) +
  geom_line(size=1.5, aes(color=ITN_Coverage, linetype=ACT_Coverage)) +
  # facet_grid(~ACT_Coverage) + 
  geom_abline() + 
  theme_minimal() + 
  labs(x="Initial", 
       y="Final",
       title="Intervention Impact, SE Asia")


