
# collate all the lookup tables used for plots that go into Gates slides

library(data.table)

main_dir <- file.path(Sys.getenv("HOME"),"Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")
out_dir <- file.path(Sys.getenv("HOME"),"Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/megatrends/gates_2019")

### "Big three" interventions, already formatted  -----------------------------------------------------
main_lut <- fread(file.path(main_dir, "lookup_full_interactions_v4.csv"))
main_lut <- unique(main_lut[, list(Site_Name, Intervention, mean_initial, mean_final)])
 
# keep one w/o ACT for DP
no_act <- main_lut[Intervention=="ITN 0.8; IRS 0.8;"]

# keep the relevant two
main_lut <- main_lut[Intervention %in% c("ITN 0.8; ACT 0.8;",
                                         "ITN 0.8; IRS 0.8; ACT 0.8;") ]

### mAB's, DHA-PQ, PEV, TBV  -----------------------------------------------------
drug_lut <- fread(file.path(main_dir, "lookup_full_interactions_experimental.csv"))
# guessing: 80% coverage, 12-month vaccines
drug_lut <- unique(drug_lut[Coverage==0.4 & !Intervention %like% "6mo", list(Site_Name, Intervention, mean_initial, mean_final)])
drug_lut[, Intervention:= gsub("_Coverage", "", Intervention)]

# keep the DP_CM version as is. For everything else,
# layer this "over" the 80/80/80 intervention from "main" 
max_current <- main_lut[Intervention== "ITN 0.8; IRS 0.8; ACT 0.8;"]

updated_drug_lut <- lapply(unique(drug_lut$Site_Name), function(site){
  this_current <- max_current[Site_Name==site]
  splObj<-smooth.spline(this_current$mean_initial, this_current$mean_final)
  this_drug <- drug_lut[Site_Name==site]
  this_drug[, combined_mean_final:=predict(splObj, mean_final)$y]
  this_drug[, combined_mean_final:=pmax(pmin(combined_mean_final, mean_final), 0)] # bound between zero and the original value
  return(this_drug)
})
updated_drug_lut <- rbindlist(updated_drug_lut)

updated_dp_lut <- lapply(unique(drug_lut$Site_Name), function(site){
  this_baseline <- no_act[Site_Name==site]
  splObj<-smooth.spline(this_baseline$mean_initial, this_baseline$mean_final)
  this_drug <- drug_lut[Site_Name==site & Intervention=="DP_CM"]
  this_drug[, combined_mean_final:=predict(splObj, mean_final)$y]
  this_drug[, combined_mean_final:=pmax(pmin(combined_mean_final, mean_final), 0)] # bound between zero and the original value
  return(this_drug)
})
updated_dp_lut <- rbindlist(updated_dp_lut)

final_drug_lut <- rbind(updated_dp_lut[Intervention=="DP_CM", list(Site_Name, Intervention, mean_initial, mean_final=combined_mean_final)],
                        updated_drug_lut[Intervention!="DP_CM", list(Site_Name, Intervention, mean_initial, mean_final=combined_mean_final)])

### ATSBs with maximized baseline ints  -----------------------------------------------------

realistic_atsbs <- fread(file.path(main_dir, "lookup_outdoor_interventions_sweep_timing.csv")) 
# keep 0.15% and 3% actual killing
realistic_atsbs <- realistic_atsbs[(Intervention=="ATSB_11percent" & Coverage==0.8) |  (Intervention=="ATSB_40percent"& Coverage==0.4)]
realistic_atsbs[, Intervention:= plyr::mapvalues(Intervention, c("ATSB_40percent", "ATSB_11percent"), c("ATSB Initial Kill:0.15%", "ATSB Initial Kill:3%"))]
realistic_atsbs[, Intervention:= paste("Baseline:80%,", Intervention)]
realistic_atsbs <- unique(realistic_atsbs[, list(Site_Name, Intervention, mean_initial, mean_final)])

elimination_atsbs <- fread(file.path(main_dir, "lookup_atsbs_to_elim.csv"))
elimination_atsbs <- elimination_atsbs[ATSB_True_Initial_Killing %in% c(0.15, 0.25)]
elimination_atsbs[, Intervention:= paste0("Baseline:80%, ", "ATSB Initial Kill:", ATSB_True_Initial_Killing*100, "%")]
elimination_atsbs <- unique(elimination_atsbs[, list(Site_Name, Intervention, mean_initial, mean_final)])

### ATSBs with a range of baseline ints  -----------------------------------------------------
atsb_runs <- c("MAP_For_Symposium_ATSB_Higher_Existing_Intervention.csv", 
               "MAP_For_Symposium_ATSB_Lower_Intervention.csv",
               "MAP_For_Symposium_ATSB_Lower_Existing_Intervention.csv",
               "MAP_For_Symposium_ATSB_No_Existing_Intervention.csv")

initial <- fread(file.path(main_dir, "../initial/MAP_II_New_Sites_Burnin.csv"))
prelim_atsb <- rbindlist(lapply(atsb_runs, function(fname){fread(file.path(main_dir, fname))}))

other_atsbs <- merge(prelim_atsb[ATSB_Initial_Effect==0.01], initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
other_atsbs[, Run_Number:=factor(Run_Number)]
other_atsbs[, Intervention:= paste0("Baseline:", ITN_Coverage*100, "%, ", "ATSB Initial Kill:", ATSB_Initial_Effect*100, "%")]
other_atsbs[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]
other_atsbs[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]
other_atsbs <- unique(other_atsbs[, list(Site_Name, Intervention, mean_initial, mean_final)])

### Combine all  -----------------------------------------------------

full_lut <- rbind(main_lut, final_drug_lut, realistic_atsbs, elimination_atsbs, other_atsbs)
full_lut <- full_lut[!Site_Name %in% c("bajonapo", "karen")]
write.csv(full_lut, file=file.path(out_dir, "full_lookup.csv"), row.names=F)
