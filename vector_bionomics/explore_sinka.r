
###############################################################################
## Explore the vector databases from Sinka et al 2016 for possible inclusion
## of bionomics (esp anthropophily, endophily) into IDM framework
################################################################################


library(data.table)
library(ggplot2)


rm(list=ls())

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/data/MAP/vector_database")

afr <- fread(file.path(main_dir, "Bionomics Africa.csv"))
afr_subset <- afr[, list(source_id, country, site, lat, long, area_type, insecticide_control, month_start, month_end, year_start, year_end,
                         species, indoor_biting_n, indoor_biting_total, indoor_biting,
                         outdoor_biting_n, outdoor_biting_total, outdoor_biting, indoor_outdoor_biting_units)]
afr_subset <- afr_subset[complete.cases(afr_subset)]


