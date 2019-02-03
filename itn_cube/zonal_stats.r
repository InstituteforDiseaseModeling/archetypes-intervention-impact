## -----------------------------------------------------------------------------------------------------------------
# ITN Cube
# zonal_aggregation.r
# 
# Amelia Bertozzi-Villa (adapted from Daniel J Weiss), Institute for Disease Modeling, University of Oxford
# Feb 2019
# 
# Zonal stats script: Use an administrative boundary surface and a population surface to aggregate rasters to admin levels.
# 
## -----------------------------------------------------------------------------------------------------------------------

library(data.table)
library(raster)

rm(list=ls())

## ------------ set up base directories and database mounts------------------
desktop.dir <- ifelse(Sys.getenv('USERPROFILE')=='', Sys.getenv('HOME'))
out.dir <- file.path(desktop.dir, 'Dropbox (IDM)/Malaria Team Folder/projects/map_itn_cube')
zed.root <- '/Volumes/map_data/'
sam.root <- '/Volumes/map_pit/sam/bld1/bras2280/ITNcube/'

## ------------ define input and output filenames and years of interest ------------------
template.filename    <- file.path(zed.root, 'GBD2017/Processing/Stages/09a_Africa_Cubes/Checkpoint_Outputs/Summary_20181122/Africa/summaries/rasters/PfPR_mean/pr_2000_rmean_Africa.tif')

raster.path          <- sam.root
# raster.path          <- file.path(zed.root, 'cubes/5km/ITN/')
# filename.string      <- '.ITN.use.yearavg.new.adj.tif'
output.filename      <- file.path(out.dir, 'test_agg_ACC.csv')
label <- "ITN Access"

pop.path.and.prefix  <- file.path(zed.root, 'GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Ouput_Pop_Unmasked_5K/ihme_corrected_frankenpop_All_Ages_3_') 
pop.suffix           <- '.tif'  # the filename characters and extension information that follows the year

zone.path            <- file.path(zed.root, 'master_geometries/Admin_Units/Global/MAP/2018/MG_5K/Rasters/')
zone.string          <- 'admin2018_0_MG_5K'
zone.filename        <- file.path(zone.path, paste0(zone.string,'.tif'))
zone.data.filename   <- file.path(zone.path, paste0(zone.string, '_config_data.csv'))
start.year           <- 2000
end.year             <- 2016
end.replicated.year  <- 2017
n.years <- end.year - start.year + 1

## ------------ set up templates, establish extents and resolutions ------------------

# function to align resolutions between two rasters
align_res <- function(rast, template.rast){
  if (!identical(res(rast), res(template.rast))) {
    rast  <- resample(rast, template.rast, method = 'ngb')
  }
  return(rast)
}

# standardize extent and resolution between template and admin filenames
admin <- raster(zone.filename)
admin.reference <- fread(zone.data.filename)
template <- raster(template.filename)

e <- extent(template)
admin <- crop(admin,e)
admin <- align_res(admin, template)

# Define the set of countries to analyze 
endemic.vec         <- unique(admin)
endemic.vec         <- endemic.vec[!endemic.vec < 1] # Remove the negative and zero values (likely no-data) 
omit.list           <- c(1993, 2029, 2042, 2046, 2085, 2086, 2114, 2135, 2160, 2162, 2169, 2174, 2185, 2192, 2204, 2225) # for Africa only - hack addition
endemic.vec         <- endemic.vec[!endemic.vec %in% omit.list]    
n.endemic.countries <- length(endemic.vec)
results.template <- admin.reference[uid %in% endemic.vec] # basis for the full dataset

if (nrow(results.template)!=n.endemic.countries){
  stop("Length of output table does not align with number of admin units in raster!")
}

# main aggregation
full.results <- lapply(start.year:end.year, function(year){
  print(year)
  input.filename <- paste(raster.path, "ITN_", year, ".USE.tif", sep='') 
  # input.filename <- paste(raster.path, year, filename.string, sep='') 
  input <- raster(input.filename)
  input <- crop(input, e)
  input <- align_res(input, template)
  
  pop.filename   <- paste(pop.path.and.prefix, year, pop.suffix, sep='')
  pop <- raster(pop.filename)
  pop   <- crop(pop, e)
  pop <- align_res(pop, template)
  
  count.raster <- input * pop
  
  # calculate zonal stats
  full.zonal <- data.table(zonal(count.raster, admin, fun='sum'))
  pop.zonal  <- data.table(zonal(pop, admin, fun='sum'))
  
  rate.zonal <- merge(full.zonal[, list(uid=zone, input_val=sum)],
                      pop.zonal[, list(uid=zone, pop=sum)],
                      by="uid", all=T)
  
  annual.results.table <- merge(results.template, 
                                rate.zonal[, list(uid, year=year, rate=input_val/pop)],
                                by="uid", all.x=T)
  
  return(annual.results.table)
})

# extend to future years 
extended.results <- lapply((end.year+1):end.replicated.year, function(new.year){
  print(new.year)
  extended.year <- copy(full.results[[n.years]])
  extended.year[, year:=new.year]
  return(extended.year)
})

# combine and reshape
full.results <- rbindlist(append(full.results, extended.results))

# remove countries with all zeros or NAs
full.results[, sum:=sum(rate), by="uid"]
full.results <- full.results[sum>0]
full.results[, sum:=NULL]

time_series <- ggplot(full.results, aes(x=year, y=rate, color=ISO3)) +
                geom_line() +
                facet_wrap(~ISO3) +
                theme(legend.position = "none") +
                labs(x="Year", y="Annual Use (Sam's Folder)")
print(time_series)

full.results.wide <- dcast(full.results, uid + ISO3 ~ year, value.var="rate")
write.csv(full.results.wide, file=output.filename, row.names=F)


