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

analysis_subdir <- "20200506_reextract_20191009_mega_era5_new_arch"
analysis_metric <- "prev"
final_day <- 365
suffix <- ""

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
out_dir <- file.path(main_dir,"results", "clean")
dir.create(out_dir, recursive = T, showWarnings = F)


# read in data
in_dir <- file.path(main_dir, "results", "raw", analysis_metric)
fnames <- list.files(in_dir)
initial <- fread(file.path(in_dir, fnames[fnames %like% "Burnin"]))
final <- fread(file.path(in_dir, fnames[fnames %like% "Int"]))

# in some cases, the final year will be duplicated, with "0" in the prevalence column. Flag and remove these.
# (should only happen for the older prevalence datasets)
remove_duplicates <- function(this_dt){
  dt_names <- names(this_dt)[!names(this_dt) %like% "prev"]
  this_dt[, count:=seq_len(.N), by=dt_names]
  this_dt <- this_dt[count==1]
  this_dt[, count:=NULL]
  return(this_dt)
}

if (analysis_metric=="prev"){
  initial <- remove_duplicates(initial)
  final <- remove_duplicates(final)
}

initial <- initial[day==max(initial$day)]
initial[, day:=NULL]

final <- final[day==final_day]
final[, day:=NULL]

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
if (analysis_metric=="prev"){
  base_metrics <- "prev"
}else if (analysis_metric=="inc"){
  base_metrics <- c("inc", "severe_inc")
}

initial_cols <- paste0("initial_", base_metrics)
final_cols <- paste0("final_", base_metrics)

colnames_to_keep <- c("Run_Number", "Site_Name", "x_Temporary_Larval_Habitat", final_cols)

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
int_impact <- merge(final, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all.x=T)

# find intervention counts and labels
int_list[, cov:=cov/100]
int_list_wide <- dcast.data.table(int_list, int_id~int, value.var="cov")
int_impact <- merge(int_impact, int_list_wide, by=unique_ints, all=T)
int_impact <- merge(int_impact, int_labels, by="int_id", all=T)
int_impact[, label:=factor(label, levels=int_labels$label)]

setcolorder(int_impact, c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat", "int_id", "label", initial_cols, final_cols, unique_ints))
int_impact <- int_impact[order(int_id, Site_Name, x_Temporary_Larval_Habitat, Run_Number)]

# test plots
int_impact[, transmission_intensity:=round(log10(x_Temporary_Larval_Habitat), 2)]


for(this_col in c(initial_cols, final_cols)){
  setnames(int_impact, this_col, "to_agg")
  int_impact[, new_col:=mean(to_agg), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)] 
  setnames(int_impact, "new_col", paste0("mean_", this_col))
  
  if (this_col %like% "final"){
    int_impact[, new_col:=min(to_agg), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)] 
    setnames(int_impact, "new_col", paste0("min_", this_col))
    int_impact[, new_col:=max(to_agg), by=list(Site_Name, x_Temporary_Larval_Habitat, int_id)] 
    setnames(int_impact, "new_col", paste0("max_", this_col))
  }
  setnames(int_impact, "to_agg", this_col)
}

# ggplot(int_impact[Site_Name==10 & transmission_intensity>0 & int_id<125 & irs==0 & Run_Number==1], aes(x=al_cm, y=itn, fill=mean_final_severe_inc)) +
#   geom_tile() +
#   # geom_text(aes(label=round(value, 2))) +
#   scale_fill_distiller(name="Incidence", palette = "Spectral") +
#   facet_wrap(.~transmission_intensity) +
#   theme_minimal() +
#   scale_x_continuous(breaks=seq(0, 0.8, 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) +
#   scale_y_continuous(breaks=seq(0, 0.8, 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) +
#   labs(x="Effective Treatment",
#        y="ITN Coverage")


int_impact <- melt(int_impact, id.vars = c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat", "int_id", "label", unique_ints))
for (this_metric in base_metrics){
  int_impact[variable %like% paste0("initial_", this_metric) | variable %like% paste0("final_", this_metric), metric:=this_metric]
}

int_impact[, variable:=gsub("(.*[initial|final])_.*", "\\1", variable)]

for_cast <- paste(paste(c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat", "int_id", "label", unique_ints, "metric"), collapse=" + "),
                  "~ variable") 

int_impact <- dcast.data.table(int_impact, formula(for_cast))


# summarize 
summary_colnames <- c("Site_Name", "int_id", "label", "x_Temporary_Larval_Habitat", "metric",
                      "initial", "final", "mean_initial", "mean_final", "min_final", "max_final", unique_ints)
summary <- unique(int_impact[, ..summary_colnames])

# save
write.csv(int_impact, file=file.path(out_dir, paste0("full_impact_", analysis_metric, suffix, ".csv")), row.names=F)
write.csv(summary, file=file.path(out_dir, paste0("summary_impact_",analysis_metric, suffix, ".csv")), row.names=F)





