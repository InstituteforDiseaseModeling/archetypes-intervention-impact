library(data.table)
library(ggplot2)
library(gridExtra)

rm(list=ls())

main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables")

initial <- fread(file.path(main_dir, "initial_prev.csv"))
#explore distribution of initial prevalences
ggplot(initial, aes(x=log10(x_Temporary_Larval_Habitat), y=initial_prev)) +
  geom_point() +
  facet_wrap(~Site_Name)

# final <- fread("final_prev.csv")
final_outbreak <- fread(file.path(main_dir, "final_prev_moine_importations.csv"))

final_outbreak <- final_outbreak[, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, ITN_Coverage,
                                        type="Outbreak Fraction", var=outbreak_fraction, final_prev)]
final_itns <- fread(file.path(main_dir, "final_prev_moine_itn_dists.csv"))
final_itns <- final_itns[, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, ITN_Coverage,
                                type="ITN Distributions", var=ITN_Distributions, final_prev)]

final <- rbind(final_outbreak, final_itns)

lookup <- merge(final, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"))
lookup[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, ITN_Coverage, type, var)]
lookup[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, ITN_Coverage, type, var)]
lookup[, ITN_Coverage:=factor(ITN_Coverage)]
lookup[, Run_Number:=factor(Run_Number)]
#write.csv(lookup, "full_lookup.csv")



outbreak <-ggplot(lookup[type=="Outbreak Fraction"], aes(x=initial_prev, y=final_prev, color=ITN_Coverage)) +
            geom_line(aes(group=interaction(ITN_Coverage, Run_Number)), alpha=0.75) + 
            geom_line(aes(x=mean_initial, y=mean_final), size=1.5) +
            theme_minimal() + 
            labs(x="Initial PfPR",
                 y="Final PfPR",
                 title="Outbreak Fractions") + 
            facet_grid(type~var)

itn_dist <- ggplot(lookup[type=="ITN Distributions"], aes(x=initial_prev, y=final_prev, color=ITN_Coverage)) +
              geom_line(aes(group=interaction(ITN_Coverage, Run_Number)), alpha=0.75) + 
              geom_line(aes(x=mean_initial, y=mean_final), size=1.5) +
              theme_minimal() + 
              labs(x="Initial PfPR",
                   y="Final PfPR",
                   title="# of ITN Distributions") + 
              facet_grid(type~var)

layout <- rbind(c(1,1,1),
                c(2,2,2))

print("saving")
full_plot <- grid.arrange(outbreak, itn_dist, layout_matrix=layout)
