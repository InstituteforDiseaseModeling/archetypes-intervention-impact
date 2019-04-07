library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
theme_set(theme_minimal(base_size = 18))

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")
get_smooth <- function(x, y){
  lo <- loess(y[y>0]~x[y>0])
  predictions <- c(y[y==0], predict(lo))
  return(pmax(predictions, rep(0, length(predictions))))
}

anthro_endo_map <- data.table(Site_Name=c("aba", "kananga", "kasama", "djibo", "gode", "moine", "bajonapo", "karen"),
                              anthro=c(74.45, 65.02, 79.04, 76.6, 75, 75.78, 50, 50),
                              endo=c(80, 85, 80.38, 55.6, 50, 52.73, 60, 24.6),
                              map_color=c("#00A08A", "#D71B5A", "#F2AD00", "#F98400", "#902E57", "#5392C2", "#7DB548", "#8971B3"))
anthro_endo_map[, human_indoor:= (anthro*endo)/100]

atsb_runs <- c("MAP_For_Symposium_ATSB_Higher_Existing_Intervention.csv", 
               "MAP_For_Symposium_ATSB_Lower_Intervention.csv",
               "MAP_For_Symposium_ATSB_Lower_Existing_Intervention.csv",
               "MAP_For_Symposium_ATSB_No_Existing_Intervention.csv")

initial <- fread(file.path(main_dir, "../initial/MAP_II_New_Sites_Burnin.csv"))
prelim_data <- rbindlist(lapply(atsb_runs, function(fname){fread(file.path(main_dir, fname))}))

all_data <- merge(prelim_data[ATSB_Initial_Effect<=0.05], initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
all_data[, Run_Number:=factor(Run_Number)]
all_data[, Intervention:= paste0("Baseline:", ITN_Coverage*100, "%, ", "ASTB Initial Kill:", ATSB_Initial_Effect*100, "%")]

all_data[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]
all_data[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]

all_data = merge(all_data, anthro_endo_map, by="Site_Name", all.x=T)

minmaxes <- all_data[, list(mean_initial=unique(mean_initial),
                             min_final=min(final_prev),
                             max_final=max(final_prev)),
                      by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]



all_data <- merge(all_data, minmaxes, by=c("Site_Name", "x_Temporary_Larval_Habitat", "Intervention", "mean_initial"), all=T)
all_data[, human_indoor:=as.factor(human_indoor)]

these_colors <- unique(all_data[!Site_Name %in% c("karen", "bajonapo"), list(human_indoor, map_color)])
these_colors <- these_colors[order(human_indoor)]$map_color

# pdf(file.path(Sys.getenv("HOME"), "Desktop", "baseline_60.pdf"), width=7, height=5)
ggplot(all_data[ITN_Coverage==0.6 & ATSB_Initial_Effect==0 & !Site_Name %in% c("karen", "bajonapo")], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.25)+
  geom_ribbon(aes(ymin=min_final, ymax=max_final, fill=human_indoor, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=human_indoor, group=Site_Name), size=2) +
  scale_color_manual(values=these_colors, name="Indoor Biting %") + 
  scale_fill_manual(values=these_colors, name="Indoor Biting %") +
  theme(legend.position="left") + 
  coord_fixed() +
  labs(x="Initial PfPR",
       y="Final PfPR",
       title="60% ITN, IRS, ACT Coverage")
# graphics.off()



pdf(file.path(Sys.getenv("HOME"), "Desktop", "baseline_60.pdf"), width=7, height=5)
ggplot(all_data[ITN_Coverage==0.6 & ATSB_Initial_Effect==0], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.25)+
  geom_ribbon(aes(ymin=min_final, ymax=max_final, fill=human_indoor, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=human_indoor, group=Site_Name), size=2) +
  scale_color_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") + 
  scale_fill_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  theme(legend.position="left") + 
  coord_fixed() +
  labs(x="Initial PfPR",
       y="Final PfPR",
       title="60% ITN, IRS, ACT Coverage")
graphics.off()

pdf(file.path(Sys.getenv("HOME"), "Desktop", "atsb_5.pdf"), width=7, height=5)
ggplot(all_data[ITN_Coverage==0 & ATSB_Initial_Effect==0.05], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.25)+
  geom_ribbon(aes(ymin=min_final, ymax=max_final, fill=human_indoor, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=human_indoor, group=Site_Name), size=2) +
  scale_color_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") + 
  scale_fill_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  theme(legend.position="left") + 
  coord_fixed() +
  labs(x="Initial PfPR",
       y="Final PfPR",
       title="ATSBs, 5% Initial Killing")
graphics.off()

pdf(file.path(Sys.getenv("HOME"), "Desktop", "baseline_0.4_atsb_5.pdf"), width=7, height=5)
ggplot(all_data[ITN_Coverage==0.4 & ATSB_Initial_Effect==0.05], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.25)+
  geom_ribbon(aes(ymin=min_final, ymax=max_final, fill=human_indoor, group=Site_Name), alpha=0.25) +
  geom_line(aes(color=human_indoor, group=Site_Name), size=2) +
  scale_color_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") + 
  scale_fill_manual(values=brewer.pal(8, "Spectral"), name="Indoor Biting %") +
  theme(legend.position="left") + 
  coord_fixed() +
  labs(x="Initial PfPR",
       y="Final PfPR",
       title="40% ITN, IRS, ACT Coverage + \n ATSBs, 5% Initial Killing")
graphics.off()


