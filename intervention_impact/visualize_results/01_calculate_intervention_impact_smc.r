###############################################################################################################
## 01_calculate_intervention_impact_SMC.r
## Amelia Bertozzi-Villa
## December 2019
## 
## 
##############################################################################################################


library(data.table)
library(ggplot2)

rm(list=ls())

analysis_subdir <- "20200423_test_smc"
main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
out_dir <- file.path(main_dir,"results", "clean")
dir.create(out_dir, recursive = T, showWarnings = F)

# read in data
in_dir <- file.path(main_dir, "results", "raw")
fnames <- list.files(in_dir)
initial <- fread(file.path(in_dir, fnames[fnames %like% "Burnin"]))
final <- fread(file.path(in_dir, fnames[fnames %like% "Int"]))
final <- final[day==365]
final[, day:=NULL]

# read original intervention specs
int_list <- fread(file.path(main_dir, "input", "interventions.csv"))
int_list[, cov:=cov/100]
unique_ints <- unique(int_list$int)
if (length(unique(int_list$start_day))==1 & unique(int_list$start_day)[1]==0){
  int_list[, start_day:=NULL]
}
setnames(int_list, c("start_day", "cov", "max_age"), c("Start", "Coverage", "Max_Age"))
int_list <- melt(int_list, id.vars=c("int_id", "int"))
int_list[, variable:=paste0(int, "_", variable)]
int_list[, variable:=gsub("itn", "ITN", variable)]
int_list[, variable:=gsub("smc", "SMC", variable)]
int_list <- dcast.data.table(int_list, int_id ~ variable, value.var="value")
# 
# # add labels -- TODO
# int_labels <- rbindlist(lapply(unique(int_list$int_id), function(this_id){
#   subset <- int_list[int_id==this_id]
#   label <- ""
#   for (int_idx in 1:nrow(subset)){
#     label <- paste0(label, subset[int_idx, int], ": ", subset[int_idx, cov], "%")
#     if(int_idx<nrow(subset)){
#       label <- paste0(label, ", ")
#     }
#   }
#   return(data.table(int_id=this_id, label=label))
# }))

# set of modifications to create a separate "coverage" column for each intervention

# collate into a single dataframe
int_impact <- merge(initial, final, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)

# find intervention counts and labels

int_impact <- merge(int_impact, int_list, by=intersect(names(int_list), names(int_impact)), all=T)
# int_impact <- merge(int_impact, int_labels, by="int_id", all=T)
# int_impact[, label:=factor(label, levels=int_labels$label)]

# setcolorder(int_impact, c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat", "int_id", "label", "initial_prev", "final_prev", unique_ints))
int_impact <- int_impact[order(int_id, Site_Name, x_Temporary_Larval_Habitat, Run_Number)]
int_impact[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]
int_impact[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]
int_impact[, min_final:= min(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]
int_impact[, max_final:=max(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]

# summarize 
# summary_colnames <- c("Site_Name", "int_id", "label", "x_Temporary_Larval_Habitat", "mean_initial", "mean_final", "min_final", "max_final", unique_ints)
# summary <- unique(int_impact[, ..summary_colnames])

# save
write.csv(int_impact, file=file.path(out_dir, "full_impact.csv"), row.names=F)
# write.csv(summary, file=file.path(out_dir, "summary_impact.csv"), row.names=F)


# plot

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

sensitivity_dir <- file.path(main_dir, "../20191218_site_sensitivity")
colormap <- fread(file.path(sensitivity_dir, "results", "clean", "cluster_color_map.csv"))
colormap <- colormap[order(ns_order)]
palette <- colormap$color

to_plot <- unique(int_impact[Site_Name<11, list(int_id, cluster=Site_Name, x_Temporary_Larval_Habitat, ITN_Coverage, SMC_Coverage, SMC_Max_Age, SMC_Start, 
                                    mean_initial, mean_final, min_final, max_final)])
to_plot <- merge(to_plot, colormap)
to_plot[, cluster_label:= factor(ns_order, labels=colormap$name)]
to_plot[, int_id:=as.factor(int_id)]

ggplot(to_plot, aes(x=mean_initial, y=mean_final, color=as.factor(SMC_Coverage), fill=cluster_label)) +
  geom_abline() + 
  # geom_ribbon(aes(ymin=min_final, ymax=max_final), alpha=0.25) + 
  geom_line(alpha=0.8, aes(linetype=as.factor(SMC_Max_Age))) +
  # scale_color_manual(values=palette) + 
  # cale_fill_manual(values=palette) + 
  facet_wrap(~cluster_label) + 
  theme_minimal() + 
  # theme(legend.position = "none") + 
  labs(x="Initial PfPR",
       y="Final PfPR",
       title="")





