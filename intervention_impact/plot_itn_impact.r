
library(data.table)
library(ggplot2)

base_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/prelim_itn_sweeps")

start_prev <- fread(file.path(base_dir, "prelim_moine_burnins.csv"))
setnames(start_prev, c("Time Of Report", "x_Temporary_Larval_Habitat", "PfPR_2to10"),
         c("report_date", "x_temp", "start_prev"))
final_prev <- fread(file.path(base_dir, "prelim_moine_itns.csv"))
setnames(final_prev, c("Time Of Report", "x_Temporary_Larval_Habitat", "add_ITN_age_season.coverage_all", "PfPR_2to10"),
         c("report_date", "x_temp", "itn_cov", "final_prev"))
final_prev <- final_prev[report_date==max(report_date)]

prev <- merge(start_prev[, list(x_temp, Run_Number, start_prev)], 
              final_prev[, list(x_temp, Run_Number, itn_cov, final_prev)],
              by=c("x_temp", "Run_Number"),
              all=T)

means <- prev[, list(start_prev=mean(start_prev), final_prev=mean(final_prev)), by=list(x_temp, itn_cov)]


png(file.path(base_dir, "prelim_itn_effect.png"), width=800, height=700, res=140)
ggplot(means, aes(x=start_prev, y=final_prev)) +
  geom_line(aes(color=factor(itn_cov)), size=1.5) +
  theme_minimal() +
  # theme(legend.position = "none") +
  labs(title="2-year effect of ITN campaign, Mozambique-like",
       x="Start Prev",
       y="End Prev")
graphics.off()

