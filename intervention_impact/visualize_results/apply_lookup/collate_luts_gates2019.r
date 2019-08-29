
# collate all the lookup tables used for plots that go into Gates slides

library(data.table)

main_dir <- file.path(Sys.getenv("HOME"),"Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")
# "big three" interventions, already formatted
 
main_lut <- fread(file.path(main_dir, "lookup_full_interactions_v4.csv"))
main_lut <- unique(main_lut[, list(Site_Name, Intervention, mean_initial, mean_final)])
 
# mAB's, DHA-PQ, PEV, TBV
drug_lut <- fread(file.path(main_dir, "lookup_full_interactions_experimental.csv"))
drug_lut <- unique(drug_lut[, list(Site_Name, Intervention, Coverage, mean_initial, mean_final)])
# todo: ask Pete which ones he used

# atsb's with max baseline int:
realistic_atsbs <- fread(file.path(main_dir, "lookup_outdoor_interventions_sweep_timing.csv")) 
# keep 0.15% and 3% actual killing
realistic_atsbs <- realistic_atsbs[(Intervention=="ATSB_11percent" & Coverage==0.8) |  (Intervention=="ATSB_40percent"& Coverage==0.4)]
#todo:  rename for clarity

elimination_atsbs <- fread(file.path(main_dir, "lookup_atsbs_to_elim.csv"))


# atsb's with a range of baseline ints:
atsb_runs <- c("MAP_For_Symposium_ATSB_Higher_Existing_Intervention.csv", 
               "MAP_For_Symposium_ATSB_Lower_Intervention.csv",
               "MAP_For_Symposium_ATSB_Lower_Existing_Intervention.csv",
               "MAP_For_Symposium_ATSB_No_Existing_Intervention.csv")

initial <- fread(file.path(main_dir, "../initial/MAP_II_New_Sites_Burnin.csv"))
prelim_data <- rbindlist(lapply(atsb_runs, function(fname){fread(file.path(main_dir, fname))}))

other_atsbs <- merge(prelim_data, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
other_atsbs[, Run_Number:=factor(Run_Number)]
other_atsbs[, Intervention:= paste0("Baseline:", ITN_Coverage*100, "%, ", "ATSB Initial Kill:", ATSB_Initial_Effect*100, "%")]
other_atsbs[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]
other_atsbs[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]




