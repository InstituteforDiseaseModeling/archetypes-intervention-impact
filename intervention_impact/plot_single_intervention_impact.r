library(data.table)
library(ggplot2)

rm(list=ls())

initial <- fread("initial_prev.csv")
final <- fread("final_prev.csv")
final <- fread("final_prev_kananga.csv")

by_vars <- c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat")
mean_vars <- c(by_vars[by_vars!="Run_Number"], "ITN_Coverage", "ITN_Distributions")

lookup <- merge(initial, final, by=by_vars)
lookup[, mean_initial:= mean(initial_prev), by=mean_vars]
lookup[, mean_final:=mean(final_prev_kananga), by=mean_vars]
lookup[, ITN_Coverage:=factor(ITN_Coverage)]
lookup[, Run_Number:=factor(Run_Number)]
#write.csv(lookup, "full_lookup.csv")

#explore distribution of initial prevalences
ggplot(initial, aes(x=log10(x_Temporary_Larval_Habitat), y=initial_prev)) + 
  geom_point() +
  facet_wrap(~Site_Name)



ggplot(lookup, aes(x=initial_prev, y=final_prev_kananga, color=ITN_Coverage)) +
    geom_line(aes(group=interaction(ITN_Coverage, Run_Number)), alpha=0.75) + 
    geom_line(aes(x=mean_initial, y=mean_final), size=1.5) +
    facet_wrap(~ITN_Distributions)