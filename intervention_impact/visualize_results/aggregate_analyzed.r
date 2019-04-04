library(data.table)
library(ggplot2)

## just a scratch plotting script at this point, really need to tighten it up


rm(list=ls())

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")

experimental_results <- T 
anthro_endo_map <- data.table(Site_Name=c("aba", "kananga", "kasama", "djibo", "gode", "moine", "bajonapo", "karen"),
                              anthro=c(74.45, 65.02, 79.04, 76.6, 75, 75.78, 50, 50),
                              endo=c(80, 85, 80.38, 55.6, 50, 52.73, 60, 24.6))
anthro_endo_map[, human_indoor:= (anthro*endo)/100]


initial <- fread(file.path(main_dir, "../initial/MAP_II_New_Sites_Burnin.csv"))
prelim_data <- fread(file.path(main_dir, "MAP_For_Symposium_ATSB_No_Existing_Intervention.csv"))

prelim_data[, Intervention:=""]


if (experimental_results==T){
  
  # prelim_data[, tbv_half_life:=TBV_Waning_Config_Decay_Time_Constant*log(2)]
  # prelim_data[, pev_half_life:=PEV_Waning_Config_Decay_Time_Constant*log(2)]
  # 
  all_data <- prelim_data[, list(Site_Name, x_Temporary_Larval_Habitat, Run_Number,
                                 # IRS_Coverage, ITN_Coverage, AL_CM_Coverage=CM_Coverage,
                                 # Start_Day=ITN_Start,
                                 # mAB = ifelse(PEV_Waning_Config_class=="WaningEffectBox", PEV_Coverage, NA),
                                 # PEV_6mo = ifelse(pev_half_life==182, PEV_Coverage, NA),
                                 # PEV_12mo = ifelse(pev_half_life==365, PEV_Coverage, NA),
                                 # TBV_6mo = ifelse(tbv_half_life==182, TBV_Coverage, NA),
                                 # TBV_12mo = ifelse(tbv_half_life==365, TBV_Coverage, NA),
                                 ATSB_0percent = ifelse(ATSB_Initial_Effect==0, ATSB_Initial_Effect, NA),
                                 ATSB_1percent = ifelse(ATSB_Initial_Effect==0.01, ATSB_Initial_Effect, NA),
                                 ATSB_3percent = ifelse(ATSB_Initial_Effect==0.03, ATSB_Initial_Effect, NA),
                                 ATSB_5percent = ifelse(ATSB_Initial_Effect==0.05, ATSB_Initial_Effect, NA),
                                 ATSB_7percent = ifelse(ATSB_Initial_Effect==0.07, ATSB_Initial_Effect, NA),
                                 ATSB_9percent = ifelse(ATSB_Initial_Effect==0.09, ATSB_Initial_Effect, NA),
                                 ATSB_11percent = ifelse(ATSB_Initial_Effect==0.11, ATSB_Initial_Effect, NA),
                                 
                                 # Larvicide = Larvicide_Coverage,
                                 # Ivermectin_7day = ifelse(Ivermectin_Duration==7, Ivermectin_Coverage, NA),
                                 # Ivermectin_14day = ifelse(Ivermectin_Duration==14, Ivermectin_Coverage, NA),
                                 # Ivermectin_30day = ifelse(Ivermectin_Duration==30, Ivermectin_Coverage, NA),
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


all_data <- merge(all_data, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)

all_data[, Run_Number:=factor(Run_Number)]

all_data[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention, Coverage)]
all_data[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention, Coverage)]

all_data[, Coverage:=factor(Coverage)]

# write.csv(all_data, file=file.path(main_dir, "lookup_outdoor_interventions_sweep_timing.csv"), row.names = F)

# to_plot <- all_data[Site_Name=="aba"]

# pdf(file.path(Sys.getenv("USERPROFILE"), 
#               "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/writing_and_presentations/ii_paper",
#               "figures/new_ints.pdf"), height=6, width=7)

all_data = merge(all_data, anthro_endo_map, by="Site_Name", all.x=T)


ggplot(all_data, aes(x=mean_initial, y=mean_final, group=Site_Name, color=human_indoor)) +
  geom_line(size=1.5) +
  geom_abline() +
  # scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
  facet_grid(. ~ Coverage) +
  theme_minimal() +
  theme(legend.position="bottom")


for (sname in unique(all_data$Site_Name)){
  # png(paste0("/Users/bertozzivill/Desktop/", sname, ".png"), height=500, width=800)
  site_plot <- ggplot(all_data[Site_Name=="moine" & Start_Day==91], aes(x=initial_prev, y=final_prev, color=Run_Number)) +
                geom_line(size=1.5) +
                geom_abline() +
                # scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
                facet_grid(Site_Name ~ Intervention) +
                theme_minimal() +
                theme(legend.position="bottom") +
                labs(title="80% 3-round Ivermectin")
  
  print(site_plot)
  # graphics.off()
}




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



