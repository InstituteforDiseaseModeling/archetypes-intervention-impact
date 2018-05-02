# Accessibility Mapping in Magude, Mozambique
# adapted from Dan Weiss' script at https://map.ox.ac.uk/research-project/accessibility_to_cities/
# 
# Amelia Bertozzi-Villa,Institute for Disease Modeling, University of Oxford
# 2018-04-05

# Citation: D.J. Weiss, A. Nelson, H.S. Gibson, W. Temperley, S. Peedell, A. Lieber, M. Hancher, E. Poyart, S. Belchior, N. Fullman, B. Mappin, U. Dalrymple, J. Rozier, 
# T.C.D. Lucas, R.E. Howes, L.S. Tusting, S.Y. Kang, E. Cameron, D. Bisanzio, K.E. Battle, S. Bhatt, and P.W. Gething. A global map of travel time to cities to assess 
# inequalities in accessibility in 2015. (2018). Nature. doi:10.1038/nature25181.
# 
rm(list=ls())

## Required Packages
require(data.table)
require(ggplot2)
library(geosphere)

# plot defaults
theme_set(theme_minimal(base_size=12))

# Root Directories
main.dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/"
grid.dir <- paste0(main.dir, "projects/Mozambique/gridded_simulation_input/")
raster.dir <- paste0(main.dir, "projects/Mozambique/incidence_calibration/rasters/")
out.dir <- paste0(raster.dir, "magude_access/")

# Input Files
hf.filename <- paste0(main.dir, "data/Mozambique/Magude/HF-GPS/Magude_HF.csv") # Just 2 columns.  Structured as [X_COORD, Y_COORD] aka [LONG, LAT].  Use a header.

# Read in the health facility names and locations
catchments <- fread(file = hf.filename)
catchments[, hf:= substr(hf, 4, 20)] # remove "CM " prefix
catchments[, hf_orig:= hf]
catchments[hf %in% c("Magude Sede", "Facazissa"), hf:="Magude-Sede-Facazissa"] 
catchments[hf %in% c("Caputine", "Panjane"), hf:="Caputine-Panjane"]
setnames(catchments, c("lng", "lat"), c("mid_x", "mid_y"))

# survey-based HF assignments:
survey.hf <- fread(paste0(out.dir, "survey_catchments.csv"))

# distance-based HF assignments
grid.points <- fread(paste0(grid.dir, "older_scripts/grid_lookup.csv"))
grid.points$catchment <- NULL

# find distances between points
distances <- distm(grid.points[,c('mid_x','mid_y')], catchments[,c('mid_x','mid_y')], fun=distVincentyEllipsoid)
colnames(distances) <- catchments$hf_orig
rownames(distances) <- grid.points$grid_cell

# assign each grid cell to its nearest hf
grid.points$Nearest <- catchments$hf[max.col(-distances)]

# keep only Chichuco, merge on survey info
chichuco <- grid.points[Nearest=="Chichuco"]
chichuco <- merge(chichuco, survey.hf[, list(grid_cell, Survey=hf_name)], all.x=T)
chichuco_distances <- distances[chichuco$grid_cell,]

# loop over factor increases, find new nearest HFs
multipliers <- seq(1, 10)

new_hfs <- lapply(multipliers, function(mult){
  new_dist <- copy(chichuco_distances)
  new_dist[,"Chichuco"] <- new_dist[,"Chichuco"]*mult

  new_df <- copy(chichuco)
  new_df$mult <- mult
  new_df$Nearest <- catchments$hf[max.col(-new_dist)]
  
  return(new_df)
})

new_hfs <- rbindlist(new_hfs)
new_hfs <- new_hfs[!is.na(Survey)]

# validation: which comes closest to the survey values?
validation <- copy(new_hfs)
validation[, correct_pred:=ifelse(Nearest==Survey, 1, 0)]
val_summary <- validation[, list(correct_percent=sum(correct_pred)/.N*100), by="mult"]
val_summary <- val_summary[order(-correct_percent)]

png(paste0(out.dir, "validation.png"))
ggplot(val_summary, aes(x=mult, y=correct_percent)) +
  geom_line(size=2) +
  labs(x="Factor",
       y="Percent Predicted Correctly",
       title="Agreement Between Distance and Survey Catchments \n For Different Multipliers on Distance-to-Chichuco")
graphics.off()

setnames(new_hfs, "Nearest", "Distance")
new_hfs <- melt(new_hfs, id.vars = c("grid_cell", "mid_x", "mid_y", "mult"), value.name="HF", variable.name="type")

png(paste0(out.dir, "distance_chichuco.png"), width=17, height=7, units="in", res=170)
ggplot(new_hfs[mult<=7], aes(x=mid_x, y=mid_y, fill=HF)) +
  geom_tile()+
  scale_fill_brewer(type="qual", palette="Set2") +
  facet_grid(type ~ mult) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.title = element_blank()) + 
  labs(title="Catchment Areas as 'Distance' to Chichuco Increases",
       x="",
       y="") 
graphics.off()











