library(data.table)
library(ggplot2)

rm(list=ls())

main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")

experimental_results <- F


initial <- fread(file.path(main_dir, "../initial/MAP_II_Sweep_Vector_Lifespan.csv"))
prelim_data <- fread(file.path(main_dir, "MAP_II_Sweep_Vector_Lifespan_Ints.csv"))

prelim_data[, Intervention:=""]

if (experimental_results==T){
  
  prelim_data[, tbv_half_life:=TBV_Waning_Config_Decay_Time_Constant*log(2)]
  prelim_data[, pev_half_life:=PEV_Waning_Config_Decay_Time_Constant*log(2)]
  
  all_data <- prelim_data[, list(Site_Name, x_Temporary_Larval_Habitat, Run_Number,
                                 DP_CM_Coverage=CM_Coverage,
                                 # DP_MDA_Coverage = MDA_Coverage,
                                 mAB_Coverage = ifelse(PEV_Waning_Config_class=="WaningEffectBox", PEV_Coverage, NA),
                                 PEV_6mo_Coverage = ifelse(pev_half_life==182, PEV_Coverage, NA),
                                 PEV_12mo_Coverage = ifelse(pev_half_life==365, PEV_Coverage, NA),
                                 TBV_6mo_Coverage = ifelse(tbv_half_life==182, TBV_Coverage, NA),
                                 TBV_12mo_Coverage = ifelse(tbv_half_life==365, TBV_Coverage, NA),
                                 final_prev)]
  
  all_data <- melt(all_data, id.vars=c("Site_Name", "x_Temporary_Larval_Habitat", "Run_Number", "final_prev"), value.name = "Coverage", variable.name = "Intervention")
  all_data <- all_data[!is.na(Coverage)]
  
  
}else{
  all_data <- copy(prelim_data)
  setnames(all_data, "CM_Coverage", "ACT_Coverage")
  
  for (int in c("ITN", "IRS", "ACT")){
    varname = paste0(int, "_Coverage")
    all_data[get(varname)!=0, Intervention:= paste0(Intervention, int, " ", get(varname), "; ") ]
  }

  
  all_data[, Intervention:= gsub(" $", "", Intervention)]
  
  all_data[Intervention=="", Intervention:="None"]
}


all_data <- merge(all_data, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat",
                                          "gambiae.Adult_Life_Expectancy"), all=T)
setnames(all_data, "gambiae.Adult_Life_Expectancy", "Vector_Life_Expectancy")
all_data <- all_data[Vector_Life_Expectancy>5]
all_data[, Run_Number:=factor(Run_Number)]
all_data[, Vector_Life_Expectancy:=factor(Vector_Life_Expectancy)]
all_data[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention, Vector_Life_Expectancy)]
all_data[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention, Vector_Life_Expectancy)]

to_plot <- all_data[Site_Name=="aba"]

# pdf(file.path(Sys.getenv("USERPROFILE"), 
#               "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/writing_and_presentations/ii_paper",
#               "figures/new_ints.pdf"), height=6, width=7)
# ggplot(to_plot, aes(x=mean_initial, y=mean_final, color=factor(Coverage))) +
#   geom_line(size=1.5) +
#   geom_abline() + 
#   scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
#   facet_wrap( ~ Intervention) +
#   theme_minimal() +
#   theme(legend.position="none")
# graphics.off()

# pdf(file.path(Sys.getenv("USERPROFILE"), 
#               "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/writing_and_presentations/ii_paper",
#               "figures/lifespan_sensitivity.pdf"), height=6, width=7)

ggplot(all_data[Site_Name=="kananga"], aes(x=mean_initial, y=mean_final, color=Vector_Life_Expectancy)) +
  geom_point(size=1.5) +
  geom_abline() + 
  # scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
  facet_wrap( ~ Intervention) +
  theme_minimal()  +
  labs(x="Mean Initial Prevalence",
       y="Mean Final Prevalence",
       title="Site: Kananga")

ggplot(all_data[Site_Name=="kananga" & Intervention=="ACT 0.8;"], aes(x=mean_initial, y=mean_final, color=Vector_Life_Expectancy)) +
  geom_point(size=1.5) +
  geom_abline() + 
  # scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
  facet_grid( ~ Vector_Life_Expectancy) +
  theme_minimal()  +
  labs(x="Mean Initial Prevalence",
       y="Mean Final Prevalence",
       title="Site: Kananga")
  # theme(legend.position="none")

ggplot(all_data[Intervention=="ACT 0.8;"], aes(x=mean_initial, color=Vector_Life_Expectancy, fill=Vector_Life_Expectancy)) +
  geom_density(alpha=0.3) +
  facet_wrap(~ Site_Name , scales="free") + 
  theme_minimal() +
  labs(x="Initial Prevalence",
       y="Final Prevalence")

# graphics.off()

write.csv(all_data, file=file.path(main_dir, "lookup_sensitivity_lifespan.csv"), row.names = F)

