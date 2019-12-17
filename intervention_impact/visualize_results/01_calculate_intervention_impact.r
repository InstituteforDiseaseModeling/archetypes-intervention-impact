###############################################################################################################
## 01_calculate_intervention_impact.r
## Amelia Bertozzi-Villa
## December 2019
## 
## 
##############################################################################################################


library(data.table)
library(ggplot2)

rm(list=ls())

analysis_subdir <- "20191008_replicate_megatrends"
main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
out_dir <- file.path(main_dir,"results", "clean")
dir.create(out_dir, recursive = T, showWarnings = F)

# read in data
initial <- fread(file.path(main_dir, "results/raw", paste0("MAP_", analysis_subdir, "_Burnin.csv")))
final <- fread(file.path(main_dir, "results/raw", paste0("MAP_", analysis_subdir, "_Intervention.csv")))

# read original intervention specs
int_list <- fread(file.path(main_dir, "input", "interventions.csv"))
unique_ints <- unique(int_list$int)
if (length(unique(int_list$start_day))==1 & unique(int_list$start_day)[1]==0){
  int_list[, start_day:=NULL]
}

# add labels
int_labels <- rbindlist(lapply(unique(int_list$int_id), function(this_id){
  subset <- int_list[int_id==this_id]
  label <- ""
  for (int_idx in 1:nrow(subset)){
    label <- paste0(label, subset[int_idx, int], ": ", subset[int_idx, cov], "%")
    if(int_idx<nrow(subset)){
      label <- paste0(label, ", ")
    }
  }
  return(data.table(int_id=this_id, label=label))
}))

# set of modifications to create a separate "coverage" column for each intervention

colnames_to_keep <- c("Run_Number", "Site_Name", "x_Temporary_Larval_Habitat", "final_prev")

# itn 
if ("itn" %in% unique_ints){
  setnames(final, "ITN_Coverage", "itn")
  colnames_to_keep <- c(colnames_to_keep, "itn")
  if ("start_day" %in% names(int_list)){
    setnames(final, "ITN_Start", "start_day")
    colnames_to_keep <- c(colnames_to_keep, "start_day")
  }
}

# irs
if ("itn" %in% unique_ints){
  setnames(final, "IRS_Coverage", "irs")
  colnames_to_keep <- c(colnames_to_keep, "irs")
}

# case management
if ("al_cm" %in% unique_ints){
  final[, al_cm:=ifelse(CM_Drug=="AL", CM_Coverage, NA)]
  colnames_to_keep <- c(colnames_to_keep, "al_cm")
}
if ("dp_cm" %in% unique_ints){
  final[, dp_cm:=ifelse(CM_Drug=="DP", CM_Coverage, NA)]
  colnames_to_keep <- c(colnames_to_keep, "dp_cm")
}

# atsbs
if ("atsb" %in% unique_ints){
  setnames(final, "ATSB_Coverage", "atsb")
  final[, atsb:=atsb*ATSB_Initial_Effect]
  colnames_to_keep <- c(colnames_to_keep, "atsb")
}

# mAbs (simulated as PEVs with box 3-month durations)
if ("mAb" %in% unique_ints){
  final[, mAb:=ifelse(PEV_Waning_Config_class=="WaningEffectBox", PEV_Coverage, NA)]
  colnames_to_keep <- c(colnames_to_keep, "mAb")
}

# pre-erythrocytic vaccines--default half-life is 12 months
if ("pev" %in% unique_ints){
  final[, pev:=ifelse(PEV_Waning_Config_class=="WaningEffectExponential", PEV_Coverage, NA)]
  colnames_to_keep <- c(colnames_to_keep, "pev")
  
  half_life <- round(unique(final[!is.na(PEV_Waning_Config_Decay_Time_Constant)]$PEV_Waning_Config_Decay_Time_Constant)*log(2))
  if (length(half_life)>1 | half_life[1]!=365){
    warning("PEV has half-lives not equal to one year!")
  }
}

# transmission-blocking vaccines--default half-life is 12 months
if ("tbv" %in% unique_ints){
  setnames(final, "TBV_Coverage", "tbv")
  colnames_to_keep <- c(colnames_to_keep, "tbv")
  
  half_life <- round(unique(final[!is.na(TBV_Waning_Config_Decay_Time_Constant)]$TBV_Waning_Config_Decay_Time_Constant)*log(2))
  if (length(half_life)>1 | half_life[1]!=365){
    warning("TBV has half-lives not equal to one year!")
  }
  
}

final <- final[, ..colnames_to_keep]

# collate into a single dataframe
int_impact <- merge(initial, final, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)

# find intervention counts and labels
int_list[, cov:=cov/100]
int_list_wide <- dcast.data.table(int_list, int_id~int, value.var="cov")
int_impact <- merge(int_impact, int_list_wide, by=unique_ints, all=T)
int_impact <- merge(int_impact, int_labels, by="int_id", all=T)
int_impact[, label:=factor(label, levels=int_labels$label)]

setcolorder(int_impact, c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat", "int_id", "label", "initial_prev", "final_prev", unique_ints))
int_impact <- int_impact[order(int_id, Site_Name, x_Temporary_Larval_Habitat, Run_Number)]
int_impact[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]
int_impact[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]
int_impact[, min_final:= min(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]
int_impact[, max_final:=max(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)]

# summarize 
summary_colnames <- c("Site_Name", "int_id", "label", "x_Temporary_Larval_Habitat", "mean_initial", "mean_final", "min_final", "max_final", unique_ints)
summary <- unique(int_impact[, ..summary_colnames])

# save
write.csv(int_impact, file=file.path(out_dir, "full_impact.csv"), row.names=F)
write.csv(summary, file=file.path(out_dir, "summary_impact.csv"), row.names=F)





