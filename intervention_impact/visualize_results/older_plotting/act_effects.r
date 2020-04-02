##################################################################################################################
## act_effects.r
## Amelia Bertozzi-Villa
## March 2020
## 
## Request from Pete: "Do you have any existing summarisations of emod that would
##                      demonstrate the marginal relationship between changing case management rates and pfpr?"

## Interpretation:
##  - how does PfPR change under different case managment rates?
##  - what is the effect of an x% increase in CM on PfPR?
##  - how does this change in different settings? At different levels of other interventions?
##################################################################################################################

library(data.table)
library(ggplot2)

rm(list=ls())

subdir <- "20191009_mega_era5_new_arch"
suffix <- "results/megatrend_ii_rasters"
main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      subdir, "results", "clean")

out_dir <- "~/Desktop"

all_impact <- fread(file.path(main_dir, "full_impact.csv"))
all_impact <- all_impact[Site_Name<11]

# subset to just interventions of interest: no new interventions
all_impact <- all_impact[int_id<=125]

# to start: al_cm *only*
impact <- all_impact[itn==0 & irs==0]
impact[, act_label:= paste0(al_cm*100, "% Effective Treatment")]

# prelim plotting
summary_impact <- unique(impact[, list(Site_Name, x_Temporary_Larval_Habitat, int_id, label, act_label, itn, irs, al_cm, mean_initial, mean_final, min_final, max_final)])


impact_curves <- ggplot(summary_impact[Site_Name==5], aes(x=mean_initial, color=act_label, fill=act_label)) +
                        geom_abline() +
                        geom_ribbon(aes(ymin=min_final, ymax=max_final), color=NA, alpha=0.5) +
                        geom_line(aes(y=mean_final), size=1) + 
                        # facet_wrap(~Site_Name) + 
                        theme_minimal() +
                        theme(legend.title = element_blank()) +
                        labs(x="Initial PfPR",
                             y="Final PfPR",
                             title="Impact Curves for Site #5")


# marginal effect: for a given change in treatment rate, what is the change in PfPR?
impact[, prev_change:=final_prev-initial_prev]

impact[, prev_change_perc:=prev_change*100]
impact[, al_cm_perc:=al_cm*100]

change_for_regression <- impact[initial_prev>0]

impact_reg <- lm(prev_change_perc ~ al_cm_perc, data=change_for_regression)

ggplot(change_for_regression, aes(x=al_cm_perc, y=prev_change_perc)) +
  geom_jitter() + 
  geom_line(aes(y=predict(impact_reg)), color="blue", size=2) + 
  theme_minimal() +
  # facet_wrap(~Site_Name) + 
  theme(legend.position="none") +
  labs(x="Change in Effective Treatment",
       y="Change in PfPR",
       title="Impact of Treatment on PfPR")




