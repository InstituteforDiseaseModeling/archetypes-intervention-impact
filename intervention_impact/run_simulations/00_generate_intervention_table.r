###############################################################################################################
## 00_generate_intervention_table.r
## Amelia Bertozzi-Villa
## October 2019
## 
## Write data table of all intervention scenarios for EMOD to run. 
##############################################################################################################

library(data.table)

rm(list=ls())

main_dir <- "/Users/bertozzivill/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact/20191008_replicate_megatrends"

##  combinatoric interventions: run every combination of intervention-coverage level
comb_int_levels <- seq(0, 80, 20)
combined <- data.table(expand.grid(itn=comb_int_levels,
                                   irs=comb_int_levels,
                                   al_cm=comb_int_levels))
combined[, int_id:=1:nrow(combined)]
all_ints <- melt(combined, id.var="int_id", value.name="cov", variable.name="int")

##  sole interventions: run individually, perhaps with some baseline

# first with an 80-80-80 baseline
baseline <- data.table(int=c("itn", "irs", "al_cm"),
                       cov=rep(80,3))

sole_ints <- data.table(expand.grid(int=c("mAb", "tbv", "pev"),
                        cov=c(40, 80)
                        )
                        )

sole_ints[, int_id := (1:nrow(sole_ints)) + max(all_ints$int_id)]

for_sole <- do.call("rbind", replicate(nrow(sole_ints), baseline, simplify = FALSE))
for_sole[, int_id:=rep(sole_ints$int_id, each=nrow(baseline))]
sole_ints <- rbind(sole_ints, for_sole)

all_ints <- rbind(all_ints, sole_ints)

# also add cm with dp
all_ints <- rbind(all_ints, data.table(int_id=max(all_ints$int_id)+1,
                                        int=c("itn", "irs", "dp_cm"),
                                        cov=rep(80, 3)))

# add atsbs with a number of different baselines
atsb_covs <- c(0.15, 3, 15, 25)
big_three_covs <- comb_int_levels

atsbs <- data.table(expand.grid(atsb=atsb_covs, cov=big_three_covs))
atsbs[, int_id := 1:nrow(atsbs)]
for_atsb <- data.table(expand.grid(int=c("itn", "irs", "al_cm"),
                                   cov=big_three_covs))
for_atsb <- dcast.data.table(for_atsb, cov~int)
atsbs <- merge(atsbs, for_atsb, by="cov")
atsbs <- melt(atsbs, id.vars = "int_id", measure.vars = c("atsb", "itn", "irs", "al_cm"),
              variable.name = "int", value.name = "cov")
atsbs[, int_id:=int_id + max(all_ints$int_id)]

all_ints <- rbind(all_ints, atsbs)

all_ints <- all_ints[order(int_id, cov)]

all_ints <- all_ints[, list(int_id, int, start_day=0, cov)]


### save
write.csv(all_ints, file=file.path(main_dir, "interventions.csv"), row.names = F)





