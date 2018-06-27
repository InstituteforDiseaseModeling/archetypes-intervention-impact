library(data.table)
library(ggplot2)

rm(list=ls())

main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/prelim_itn_sweeps")

sites <- c("karen", "moine", "bajonapo")

# lookup_all <- lapply(sites, function(site){
#   lookup <- fread(file.path(main_dir, paste0("lookup_table_", site, "_multi_int.csv")))
#   lookup[, site:=site]
#   return(lookup)
# })
# 
# lookup_all <- rbindlist(lookup_all)
# write.csv(lookup_all, file=file.path(main_dir, "lookup_multi_site.csv"), row.names=F)
# 
# lookup_all[, ITN_Coverage:=factor(ITN_Coverage)]
# lookup_all[, IRS_Coverage:=paste0("IRS: ", IRS_Coverage)]
# lookup_all[, ACT_Coverage:= factor(ACT_Coverage)]
# 
# ggplot(lookup_all, aes(x=initial, y=final)) +
#   geom_line(size=1.5, aes(color=ITN_Coverage, linetype=ACT_Coverage)) +
#   facet_grid(site~IRS_Coverage) + 
#   geom_abline() + 
#   theme_minimal() + 
#   labs(x="Initial", 
#        y="Final",
#        title="Intervention Impact")


data_all <- fread(file.path(main_dir, "full_endo.csv"))

# split into initial/final, keep only itn-only interventions
initial <- data_all[run_type %like% "initial"]

final <- data_all[run_type %like% "final"]
final <- final[ITN_Coverage==0.8 & IRS_Coverage==0 & ACT_Coverage==0]
final[, year:= year + max(initial$year)]

data_all <- rbind(initial, final)
data_all[, endo_perc:= substr(run_type, nchar(run_type)-1, nchar(run_type))]
data_all <- data_all[, list(Run_Number, x_Temporary_Larval_Habitat, endo_perc, year, PfPR_2to10)]

ggplot(data_all[Run_Number==0 & year>45], aes(x=year, y=PfPR_2to10)) +
  geom_line(aes(color=endo_perc), size=1) +
  facet_wrap(~x_Temporary_Larval_Habitat)


data_by_age <- fread(file.path(main_dir, "full_karen.csv"))
data_by_age <- data_by_age[year<53]
data_by_age <- data_by_age[ run_type=="initial" | (ITN_Coverage==0.8 & IRS_Coverage==0 & ACT_Coverage==0), 
                    list(Run_Number, x_Temporary_Larval_Habitat, year, age, Average_Population_by_Age_Bin, PfPR_by_Age_Bin)]

ggplot(data_by_age[Run_Number==0 & x_Temporary_Larval_Habitat==1000 & year>45], aes(x=year, y=PfPR_by_Age_Bin)) +
  geom_line() +
  facet_wrap(~age)


# calculate all_age prevalence by year

data_by_age[, infected:=PfPR_by_Age_Bin*Average_Population_by_Age_Bin]
total_prev <- data_by_age[, list(pop=sum(Average_Population_by_Age_Bin), infected=sum(infected)), 
                          by=list(Run_Number, x_Temporary_Larval_Habitat, year)]
total_prev[, prev:=infected/pop]


# age_dist <- data_all[year %in% c(0, 10, 20, 30, 40, 50), list(pop=mean(Average_Population_by_Age_Bin)), by=list(age, year)]
# 
# sum(age_dist[age<10 & age>=2]$pop)/sum(age_dist$pop)
# 
# 
# ggplot(data_all, aes(x=year, y=Average_Population_by_Age_Bin)) +
#   geom_line(aes(color=factor(x_Temporary_Larval_Habitat))) +
#   facet_wrap(~age, scales="free_y") +
#   theme(legend.position = "none")



# png(file.path(main_dir, "anthro_sweep.png"), width=3, height=2, units="in", res=100)
# ggplot(lookup_karen, aes(x=initial, y=final)) +
#   geom_line(size=1.5) +
#   # facet_grid(.~type) +
#   theme_minimal() +
#   labs(title="SE Asia ITN+ACT Impact",
#        x="Initial Prevalence",
#        y="Final Prevalence")
# graphics.off()



