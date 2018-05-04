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
require(gdistance)
require(abind)
require(data.table)
require(ggplot2)
require(RANN)
require(MapSuite)
require(scales)

# plot defaults
theme_set(theme_minimal(base_size=14))

# Set to TRUE if a transition matrix (T.GC) has already been created and saved.   
transition.matrix.exists.flag <- T
plot.raster <- F

# Root Directories
main.dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/"
grid.dir <- paste0(main.dir, "projects/Mozambique/gridded_simulation_input/")
raster.dir <- paste0(main.dir, "projects/Mozambique/incidence_calibration/rasters/")
out.dir <- paste0(raster.dir, "magude_access/")

if (!dir.exists(out.dir)){
  dir.create(out.dir)
}

# Input Files
friction.surface.filename <- paste0(raster.dir,"friction_surface/friction_surface_all.year.2015.tif")
hf.filename <- paste0(main.dir, "data/Mozambique/Magude/HF-GPS/Magude_HF.csv") # Just 2 columns.  Structured as [X_COORD, Y_COORD] aka [LONG, LAT].  Use a header.

# Output Files
T.filename <- paste0(out.dir, 'study.area.T.rds')
T.GC.filename <- paste0(out.dir, 'study.area.T.GC.rds')
output.raster.filename <- paste0(out.dir, 'study.area.accessibility.tif')

# Read in the health facility names and locations
catchments <- fread(file = hf.filename)
catchments[, hf:= substr(hf, 4, 20)] # remove "CM " prefix
catchments[, hf_orig:= hf]
catchments[hf %in% c("Magude Sede", "Facazissa"), hf:="Magude-Sede-Facazissa"] 
catchments[hf %in% c("Caputine", "Panjane"), hf:="Caputine-Panjane"]


## Distance: 
# read in gridded input file (with distance-based catchment areas)
orig.grid.points <- fread(paste0(grid.dir, "older_scripts/grid_lookup.csv"))
orig.grid.points[catchment %in% c("Magude-Sede", "Facazissa"), catchment:="Magude-Sede-Facazissa"]
orig.grid.points[catchment %in% c("Caputine", "Panjane"), catchment:="Caputine-Panjane"]
grid.points <- copy(orig.grid.points)
setnames(grid.points, "catchment", "Distance")


# test facebook
fb.surface.filename <- paste0(raster.dir,"facebook_pop/facebook_pop_all.tif")
fb_pop = raster(fb.surface.filename)
agg_fb = aggregate(fb_pop, fact=33, fun=sum)

fb_pop_df = copy(orig.grid.points)
grid.indices <- cellFromXY(agg_fb, grid.points[, list(mid_x, mid_y)])
fb_pop_df$fb_pop <- agg_fb[grid.indices]
write.csv(fb_pop_df, file= paste0(raster.dir, "../transmission_population_comparison/fb_pop.csv"), row.names=F)



## Survey:
# read in household health facility preferences, assign grid points to preferred hf's
hhs <- fread(paste0(main.dir, "data/Mozambique/Magude/HHs/hh_by_health_facility.csv"))
hhs <- hhs[!hf_name %in% c("CHW", "Traditional healer", "Other")]
hhs[hf_name=="Chicuco", hf_name:="Chichuco"] # for consistency with other data sources
hhs[hf_name %in% c("Magude-Sede", "Facazissa"), hf_name:="Magude-Sede-Facazissa"]
hhs[hf_name %in% c("Caputine", "Panjane"), hf_name:="Caputine-Panjane"]
hhs$House_ID <- 1:nrow(hhs)

# nearest neighbor search to assign houses to grid cells
hhs$grid_cell <- nn2(grid.points[, list(mid_x, mid_y)], hhs[, list(lng_r2, lat_r2)], k=1)$nn.idx 
hh.catch <- hhs[, list(grid_hf_count=.N), by=list(grid_cell, hf_name)]
hh.catch <- hh.catch[order(grid_cell, -grid_hf_count)]
# after sorting, the first HF per grid cell is the most-utilized. Keep only this.
hh.catch <- hh.catch[, .SD[1], by=list(grid_cell)] 

write.csv(hh.catch, file=paste0(out.dir, "survey_catchments.csv"), row.names=F)

# note: household survey apparently failed to capture a large swath of the southeast
#png(paste0(out.dir, "survey_to_grid.png"))
ggplot(orig.grid.points, aes(x=mid_x, y=mid_y)) +
  geom_tile(fill="darkturquoise") +
  geom_point(data=hhs, aes(x=lng_r2, y=lat_r2)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=2)) +
  labs(title="Simulation Grid Cells vs. Household Locations, \n Census Round 2",
       x="",
       y="")
#graphics.off()


## Friction:
#  Define the spatial template
fs1 <- raster(friction.surface.filename)

# Make the graph and the geocorrected version of the graph (or read in the latter).
if (transition.matrix.exists.flag) {
  # Read in the transition matrix object if it has been pre-computed
  T.GC <- readRDS(T.GC.filename)
} else {
  # Make and geocorrect the transition matrix (i.e., the graph)
  T <- transition(fs1, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
  saveRDS(T, T.filename)
  T.GC <- geoCorrection(T)
  saveRDS(T.GC, T.GC.filename)
}

# Loop through catchments, extract travel times to each grid cell
friction <- grid.points[, list(grid_cell, mid_x, mid_y)]

for (facility in unique(catchments$hf)){
  print(facility)
  this.point <- as.matrix(catchments[hf==facility, list(X_COORD=lng, Y_COORD=lat)])

  # Run the accumulated cost algorithm to make the accessibility map
  temp.raster <- accCost(T.GC, this.point)
  grid.indices <- cellFromXY(temp.raster, friction[, list(mid_x, mid_y)])
  friction[[facility]] <- temp.raster[grid.indices]

  # convert raster to data table for plotting
  raster.dt <- as.data.table(as.data.frame(temp.raster, xy=T))
  raster.dt$id <- rownames(raster.dt)
  raster.dt$hf <- facility

  # plot full raster
  if (plot.raster){
    RasterMap(raster.dt, id="id", xcol="x", ycol="y", variable="layer", map_colors = wpal("sky"),
              font_family = "sans", map_title=paste("Accessibility:", facility),
              map_title_font_size=14, legend_title="Minutes of Travel",
              pdf_path = paste0(out.dir, "accessibility_", facility, ".pdf"))
  }
}


# find minimum travel times
friction <- melt(friction[!is.na(Moine)], id.vars=c("grid_cell", "mid_x", "mid_y"),
                 variable.name = "hf", value.name="time_to_hf")

friction[, min_time:=min(time_to_hf), by=list(grid_cell)]
friction.mins <- friction[time_to_hf==min_time]


## Merge all catchments; reshape and plot
grid.points <- merge(grid.points, friction.mins[, list(grid_cell, Friction=hf)], by="grid_cell")
grid.points <- merge(grid.points, hh.catch[, list(grid_cell, Survey=hf_name)], by="grid_cell")

grid.points <- melt(grid.points, id.vars = c("grid_cell", "mid_x", "mid_y"),
                    variable.name = "type", value.name="catchment")

png(paste0(out.dir, "catchment_compare.png"), width=1250, height=550)
ggplot(grid.points, aes(x=mid_x, y=mid_y)) +
  geom_tile(aes(fill=catchment))+
  facet_grid(~type) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.title = element_blank()) +
  labs(title="Catchment Areas Determined by Distance to HF, \n Health-Seeking Survey, or Friction Surface",
       x="",
       y="") 
graphics.off()


## New Gridded Input Files
# generate new health-seeking file based on friction
youth.clin.max = 0.9
adult.clin.max = 0.6
severe.min = 0.8

distance.healthseek <- fread(paste0(grid.dir, "grid_all_healthseek_events.csv"))

grid.healthseek <- friction.mins[, list(grid_cell,
                                        cov_newclin_youth=rescale(min_time, c(youth.clin.max,0)),
                                        cov_newclin_adult=rescale(min_time, c(adult.clin.max,0)))]
grid.healthseek[,cov_severe_youth:= pmax(severe.min, cov_newclin_youth)]
grid.healthseek[,cov_severe_adult:= pmax(severe.min, cov_newclin_youth)]

# fill missing rows with distance.healthseek values
grid.healthseek <- merge(grid.healthseek, distance.healthseek[, list(grid_cell, event, fulldate, duration)], by="grid_cell")
grid.healthseek <- rbind(grid.healthseek, distance.healthseek[!grid_cell %in% unique(grid.healthseek$grid_cell)])
grid.healthseek <- grid.healthseek[order(grid_cell)]
setcolorder(grid.healthseek, c("grid_cell", "event", "fulldate", "duration",
                               "cov_newclin_youth", "cov_newclin_adult", 
                               "cov_severe_youth", "cov_severe_adult"))

# generate new lookup file based on friction
grid.lookup <- merge(orig.grid.points, friction.mins[, list(grid_cell, Friction=hf)], by="grid_cell", all=T)
grid.lookup <- grid.lookup[, list(grid_cell, mid_x, mid_y, catchment=ifelse(is.na(Friction), catchment, as.character(Friction)))]


# save new gridded inputs
write.csv(grid.healthseek, file=paste0(grid.dir, "grid_all_healthseek_events_friction.csv"), row.names = F)
write.csv(grid.lookup, file=paste0(grid.dir, "grid_lookup_friction.csv"), row.names=F)


