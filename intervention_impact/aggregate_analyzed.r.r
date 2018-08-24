library(data.table)
library(ggplot2)

files <- list.files(file.path(main_dir, "interactions", "version_2"), full.names = T)
all_data <- lapply(files, fread)
all_data <- rbindlist(all_data, fill=T)
all_data[, Intervention:=""]


for (int in c("ITN", "IRS", "ACT")){
  varname = paste0(int, "_Coverage")
  all_data[get(varname)!=0, Intervention:= paste0(Intervention, int, " ", get(varname), "; ") ]
}

all_data[Intervention=="", Intervention:="None"]
all_data <- merge(all_data, initial, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)

all_data[, Run_Number:=factor(Run_Number)]
all_data[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]
all_data[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, Intervention)]

write.csv(all_data, file=file.path(main_dir, "interactions", "version_2", "lookup_full_interactions_v2.csv"), row.names = F)

