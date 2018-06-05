
library(data.table)
library(ggplot2)

base_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/prelim_itn_sweeps")

start_prev <- fread(file.path(base_dir, "prelim_moine_burnins.csv"))
setnames(start_prev, c("Time Of Report", "x_Temporary_Larval_Habitat", "PfPR_2to10"),
         c("report_date", "x_temp", "start_prev"))

final_prev <- fread(file.path(base_dir, "moine_itn_discard.csv"))

# static names
setnames(final_prev, c("Time Of Report", "x_Temporary_Larval_Habitat",  "PfPR_2to10"),
         c("report_date", "x_temp", "final_prev"))

# possibly variable intervention names-- add dummy columns if needed
int_rename <- c(itn_cov = "add_ITN_age_season.coverage_all", cm_cov="Healthseek_Coverage", irs_cov="IRS_Coverage")
for (val in names(int_rename)){
  if (int_rename[[val]] %in% names(final_prev)){
    setnames(final_prev, int_rename[[val]], val)
  }else{
    final_prev[[val]] <- 0
  }
}


prev <- merge(start_prev[, list(x_temp, Run_Number, start_prev)], 
              final_prev[, list(x_temp, Run_Number, itn_cov, cm_cov, irs_cov, final_prev)],
              by=c("x_temp", "Run_Number"),
              all=T)

means <- prev[, list(start_prev=mean(start_prev), final_prev=mean(final_prev)), 
              by=list(x_temp, itn_cov, cm_cov, irs_cov)]
means[, cm_cov:= paste("CM:", cm_cov)]
means[, itn_cov:=paste("ITN:", itn_cov)]
means[, irs_cov:=paste("IRS:", irs_cov)]

#png(file.path(base_dir, "cm_itn_effect.png"), width=1000, height=700, res=140)
ggplot(means, aes(x=start_prev, y=final_prev)) +
  geom_line(aes(color=factor(itn_cov)), size=1.5) +
  theme_minimal() +
  facet_wrap(~cm_cov) +
  theme(legend.title = element_blank()) +
  labs(title="2-year effect of ITNs, Mozambique-like",
       x="Start Prev",
       y="End Prev")
# graphics.off()

