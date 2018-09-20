library(data.table)
library(ggplot2)
library(gridExtra)
library(viridis)

rm(list=ls())
theme_set(theme_minimal())

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

jaline_colors <- list(purple="#AF81AD",
                      red="#E9806C",
                      orange="#F1B657",
                      green="#B1D066",
                      teal="#89CBBF",
                      blue="#7EACD9",
                      indigo="#BAB6CE"
)

main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables")

initial <- fread(file.path(main_dir, "/initial/initial_burnin_3.csv"))
final <- fread(file.path(main_dir, "interactions", "ucsf_hs_sweep.csv"))

all_data <- merge(final, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)

all_data[, Run_Number:=factor(Run_Number)]
all_data[, Health_Seeking_Time:=factor(ACT_Daily_Prob, labels=c("Two Weeks", "One Week", "Two Days"))]
all_data[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, ACT_Daily_Prob)]
all_data[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, ACT_Daily_Prob)]

means <- unique(all_data[, list(Site_Name, x_Temporary_Larval_Habitat, mean_initial, mean_final, Health_Seeking_Time)])

site_colors <- unname(unlist(jaline_colors[c("red", "orange", "green")]))

ggplot(means, aes(x=mean_initial, y=mean_final, color=Health_Seeking_Time)) +
  geom_line(size=1.25) + 
  facet_wrap(~Site_Name) +
  scale_color_manual(values=site_colors) + 
  labs(x="Initial PfPR",
       y="Final PfPR")

write.csv(means, file=file.path(main_dir, "interactions", "lookup_ucsf_hs_sweep.csv"))
