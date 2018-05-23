# Accessibility Mapping in Northcoast, Peru
# adapted from Dan Weiss' script at https://map.ox.ac.uk/research-project/accessibility_to_cities/
# 
# Amelia Bertozzi-Villa,Institute for Disease Modeling, University of Oxford
# 2018-04-17

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

# Root Directories
main.dir <- "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/"
raster.dir <- paste0(main.dir, "projects/peru_northcoast/friction_surface/rasters/")
out.dir <- paste0(raster.dir, "accessibility/")

if (!dir.exists(out.dir)){
  dir.create(out.dir)
}

# Input Files
friction.surface.filename <- paste0(raster.dir,"friction/friction_all.year.2015.tif")
hf.filename <- paste0(main.dir, "data/Peru/AntonioQuispe_VisitingScholar2018/NorthCoast_Peru/northcoast1.csv") # Just 2 columns.  Structured as [X_COORD, Y_COORD] aka [LONG, LAT].  Use a header.

# Output Files
T.filename <- paste0(out.dir, 'study.area.T.rds')
T.GC.filename <- paste0(out.dir, 'study.area.T.GC.rds')
output.raster.filename <- paste0(out.dir, 'study.area.accessibility.tif')

# Read in the health facility names and locations
hfs <- fread(file = hf.filename)
# hfs <- hfs[year==2011 & week==1] # get unique locations

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

ecuador_border <- data.table(longitudw=-80.24389, latituds= -3.48149)
points <- as.matrix(ecuador_border[, list(X_COORD=longitudw, Y_COORD=latituds)])

# Run the accumulated cost algorithm to make the accessibility map
temp.raster <- accCost(T.GC, points)

grid.indices <- cellFromXY(temp.raster, hfs[, list(longitudw, latituds)])
hfs$minutes_to_border <- temp.raster[grid.indices]
write.csv(hfs, file=paste0(out.dir, "weely_northcoast_with_accessibility.csv"), row.names=F)

# convert raster to data table for plotting
raster.dt <- as.data.table(as.data.frame(temp.raster, xy=T))
raster.dt$id <- rownames(raster.dt)

png(paste0(out.dir, "accessibility_ecuador.png"), width=700, height=550)
ggplot(raster.dt, aes(x=x, y=y)) +
  geom_tile(aes(fill=layer)) +
  scale_fill_gradientn(colors=wpal("sky")) +
  geom_point(data=ecuador_border, aes(x=longitudw, y=latituds), size=3) +
  theme(legend.title = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title="Minutes of Travel Time to Ecuadorian Border",
       x="Longitude",
       y="Latitude")
graphics.off()

writeRaster(temp.raster, output.raster.filename)





