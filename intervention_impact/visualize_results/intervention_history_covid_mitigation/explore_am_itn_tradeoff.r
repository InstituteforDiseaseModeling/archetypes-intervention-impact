###############################################################################################################
## explore_am_itn_tradeoff.r
## Amelia Bertozzi-Villa
## April 2020
## 
## See how high other interventions need to be to make up for lowered antimalarial coverage
##############################################################################################################

library(data.table)
library(ggplot2)
library(pracma)

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

# to start: no irs
impact <- all_impact[irs==0]
impact[, al_label:= paste0(al_cm*100, "% Effective Treatment")]
impact[, itn_label:= paste0(itn*100, "% ITN Coverage")]
impact[, transmission_intensity:=round(log10(x_Temporary_Larval_Habitat), 2)]

summary_impact <- unique(impact[, list(Site_Name, transmission_intensity, int_id, itn, irs, al_cm, prev=mean_final)])


    x <- linspace(-1, 1, 11)
    y <- linspace(-1, 1, 11)
    mgrid <- meshgrid(x, y)
    Z <- mgrid$X^2 + mgrid$Y^2
    xp <- yp <- linspace(-1, 1, 101)

    method <- "linear"
    zp <- interp2(x, y, Z, xp, yp, method)
    plot(xp, zp, type = "l", col = "blue")


unique_itn <- sort(unique(summary_impact$itn))
unique_act <- sort(unique(summary_impact$al_cm))
this_site <- 8
this_orig_itn_cov <- 0.4
this_orig_act_cov <- 0.6
this_new_act_cov <- 0.48
this_prev <- 0.32

this_impact <- summary_impact[Site_Name==this_site]

impact_matrices <- lapply(unique(this_impact$transmission_intensity), function(this_intensity){
  for_interp <- this_impact[transmission_intensity==this_intensity]
  # row is itn, column is act
  for_interp <- as.matrix(dcast.data.table(for_interp, itn~al_cm, value.var = "prev")[, 2:(length(unique_itn)+1)])
  rownames(for_interp) <- unique_itn
  return(for_interp)
})
names(impact_matrices) <- unique(this_impact$transmission_intensity)


estimated_prevs <- sapply(impact_matrices, function(this_matrix){
  return(interp2(unique_itn, unique_act, this_matrix, this_orig_act_cov, this_orig_itn_cov))
}) 
prev_diffs <- abs(estimated_prevs-this_prev)
intensity_id <- names(which(prev_diffs==min(prev_diffs)))

intensity_matrix <- impact_matrices[[intensity_id]]
possible_prevs <- interp2(unique_itn, unique_act, intensity_matrix, rep(this_new_act_cov, length(unique_itn)), unique_itn)
names(possible_prevs) <- unique_itn

this_new_itn_coverage <- interp1(rev(possible_prevs), rev(unique_itn), this_prev)



# pdf("~/Desktop/itn_am_heatmaps.pdf", width=14, height=8)

for (this_site in unique(impact$Site_Name)){
  # this_heatmap <- 
    ggplot(summary_impact[Site_Name==8 & transmission_intensity>0], aes(x=al_cm, y=itn, fill=mean_final)) + 
    geom_tile() +
    geom_text(aes(label=round(mean_final, 2))) +
    scale_fill_distiller(name="Prevalence", palette = "Spectral") +
    facet_wrap(.~transmission_intensity) + 
    theme_minimal() +
    scale_x_continuous(breaks=seq(0, 0.8, 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) + 
    scale_y_continuous(breaks=seq(0, 0.8, 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) + 
    labs(x="Effective Treatment",
         y="ITN Coverage")
  # print(this_heatmap)
}

# graphics.off()


## a new take: what if we explore the space of the outputs MAP has already produced instead? 
map_results_fname <- "~/Dropbox (IDM)/Malaria Team Folder/projects/map_itn_cube/covid_intervention_scenarios/combined_counterfactual_with_deathsNATIONAL_ITN6_AM3_rakedWMR2019_FORSHP.csv"
map_results_raw <- fread(map_results_fname)
map_results_raw[, c("ISO", "Admn_Level", "ID_0", "Code_0", "Type_0",
                    "Source", "Name"):=NULL]
map_results <- melt(map_results_raw, id.vars=c("Name_0", "ISO3", "PopSUM"), variable.name = "full_var")
map_results <- map_results[!is.na(value), list(Country=Name_0,
                                               iso3=ISO3,
                                               pop=PopSUM,
                                               variable=gsub("(In|De)_.*", "\\1", full_var),
                                               statistic=gsub(".*_(.*)_.*", "\\1", full_var),
                                               scenario=gsub(".*_.*_(.*)", "\\1", full_var),
                                               value
                                               )]
map_subset <- map_results[statistic=="Ex"  &  variable=="De" & iso3=="MOZ"]

ggplot(map_subset) +
  geom_bar(data=map_subset, aes(x=scenario, y=value), stat="identity") +
  # geom_hline(data=map_subset[scenario %like% "BO"], aes(yintercept=value, color=scenario)) + 
  facet_wrap(~iso3) +
  theme_minimal()





