library(data.table)
library(ggplot2)

## response to reviewer #4 on burn-ins
theme_set(theme_minimal(base_size = 18))

rm(list=ls())

get_smooth <- function(x, y){
  
  if (max(y)<0.05){
    return(y)
  }else{
    lo <- loess(y[y>0]~x[y>0])
    predictions <- c(y[y==0], predict(lo))
    return(pmax(predictions, rep(0, length(predictions))))
  }
}

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")

new_names <- data.table(Site_Name=c("aba", "bajonapo", "djibo", "gode", "kananga", "karen", "kasama", "moine"),
                        Print_Name=c("Nigeria (teal)", "Americas", "Burkina Faso (orange)",
                                     "Ethiopia (purple)", "DRC (magenta)", "SE Asia", "Zambia (blue)",
                                     "Mozambique (yellow)"))

# long burnin
long_initial <- fread(file.path(main_dir, "../initial/MAP_Longer_Oldclimate_Burnin.csv"))
long_final <- fread(file.path(main_dir, "MAP_CM_Caitlin_Longer_Burnin_Intervention.csv"))
long_all <- merge(long_final, long_initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
long_all[, burnin_type:="40-Year Burn-in"]

short_initial <- fread(file.path(main_dir, "../initial/MAP_II_New_Sites_Burnin.csv"))
short_final <- fread(file.path(main_dir, "MAP_II_New_Sites.csv"))
short_final <- short_final[IRS_Coverage==0 & ITN_Coverage==0, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, CM_Drug, CM_Coverage, final_prev)]
short_all <- merge(short_final, short_initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
short_all[, burnin_type:="15-Year Burn-in"]

all_data <- rbind(long_all, short_all)

all_data[, Run_Number:=factor(Run_Number)]
all_data[, Coverage:=factor(CM_Coverage)]
all_data[, mean_initial:= mean(initial_prev), by=list(burnin_type, Site_Name, x_Temporary_Larval_Habitat, Coverage)]
all_data[, mean_final:=mean(final_prev), by=list(burnin_type, Site_Name, x_Temporary_Larval_Habitat, Coverage)]

minmaxes <- all_data[, list(mean_initial=unique(mean_initial),
                            mean_final=unique(mean_final),
                            min_final=min(final_prev),
                            max_final=max(final_prev)),
                     by=list(burnin_type, Site_Name, x_Temporary_Larval_Habitat, Coverage)]

minmaxes_smooth <- lapply(unique(minmaxes$Site_Name), function(site_name){
  sub_list <- lapply(unique(minmaxes$burnin_type), function(burn_name){
    sub_sub_list <- lapply(unique(minmaxes$Coverage), function(cov){
      subset <- minmaxes[Site_Name==site_name & burnin_type==burn_name & Coverage==cov]
      subset[, smooth_min:= get_smooth(mean_initial, min_final)]
      subset[, smooth_max:= get_smooth(mean_initial, max_final)]
      subset[, smooth_mean:= get_smooth(mean_initial, mean_final)]
    })
    sub_sub_list <- rbindlist(sub_sub_list)
  })
  sub_list <- rbindlist(sub_list)
})
minmaxes_smooth <- rbindlist(minmaxes_smooth)

all_data <- merge(all_data, minmaxes_smooth, by=c("burnin_type", "Site_Name", "x_Temporary_Larval_Habitat", "Coverage",
                                                  "mean_initial", "mean_final"), all=T)


all_data <- merge(all_data, new_names, by="Site_Name", all=T)
all_data[, Print_Name := factor(Print_Name, levels=c("Ethiopia (purple)",
                                                     "Burkina Faso (orange)",
                                                     "Nigeria (teal)",
                                                     "DRC (magenta)",
                                                     "Zambia (blue)",
                                                     "Mozambique (yellow)",
                                                     "Americas",
                                                     "SE Asia"
                                                     ))]

png(file.path(Sys.getenv("HOME"), "Desktop/reviewer_burnin.png"), height=500, width=600)
ggplot(all_data[Coverage==0.8], aes(x=mean_initial, y=mean_final, color=burnin_type)) +
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=burnin_type), alpha=0.25) +
  geom_line() +
  geom_abline() +
  # scale_color_manual(values=c("#E9806C", "#F1B657","#B1D066")) +
  facet_wrap(Print_Name~.) +
  theme(legend.position="bottom",
        legend.title = element_blank()) +
  labs(title="80% ACT Use Over Three Years",
       x="Initial Prevalence",
       y="Final Prevalence")
graphics.off()

