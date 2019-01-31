library(raster)

template.filename    <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\rasters\\PfPR_rmean\\PfPR_rmean_Africa_admin1_2000.tif'
#raster.path          <- 'Z:\\cubes\\5km\\ITN\\'
#filename.string      <- '.ITN.use.yearavg.new.adj.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\ITN.csv'
#raster.path          <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20180423\\Africa\\anti_malarial_access\\'
#filename.string      <- '_africa_antimalarial_access_Africa.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\ACT.csv'
#raster.path          <- 'Z:\\cubes\\5km\\IRS\\'
#filename.string      <- '.IRS.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\IRS.csv'
#raster.path          <- 'Z:\\mastergrids\\MODIS_Global\\MCD43D6_v6_BRDF_Reflectance\\EVI_v6\\5km\\Annual\\EVI_v6.'
#filename.string      <- '.Annual.mean.5km.mean.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\EVI.csv'
#raster.path          <- 'Z:\\mastergrids\\MODIS_Global\\MOD11A2_v6_LST\\LST_Night\\5km\\Annual\\LST_Night_v6.'
#filename.string      <- '.Annual.mean.5km.mean.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\LSTpm.csv'
#raster.path          <- 'Z:\\mastergrids\\MODIS_Global\\MOD11A2_v6_LST\\LST_Day\\5km\\Annual\\LST_Day_v6.'
#filename.string      <- '.Annual.mean.5km.mean.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\LSTam.csv'
#raster.path          <- 'Z:\\mastergrids\\MODIS_Global\\MOD11A2_v6_LST\\LST_DiurnalDifference\\5km\\Annual\\LST_DiurnalDiff_v6.'
#filename.string      <- '.Annual.mean.5km.mean.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\LSTdelta.csv'
#raster.path          <- 'Z:\\mastergrids\\MODIS_Global\\MCD43D6_v6_BRDF_Reflectance\\TCB_v6\\5km\\Annual\\TCB_v6.'
#filename.string      <- '.Annual.mean.5km.mean.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\TCB.csv'
#raster.path          <- 'Z:\\mastergrids\\MODIS_Global\\MCD43D6_v6_BRDF_Reflectance\\TCW_v6\\5km\\Annual\\TCW_v6.'
#filename.string      <- '.Annual.mean.5km.mean.tif' # the filename string searched for (i.e., a portion of the filename)
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\TCW.csv'
#raster.path          <- 'Z:\\mastergrids\\Other_Global_Covariates\\TemperatureSuitability\\TSI_Pf_Dynamic\\5km\\Annual\\TSI-Martens2-Pf.'
#filename.string      <- '.Annual.Mean.5km.Data.tif'
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\TSI.csv'
raster.path          <- 'Z:\\mastergrids\\Other_Global_Covariates\\Rainfall\\CHIRPS\\5km\\Annual\\chirps-v2-0.'
filename.string      <- '.Annual.sum.5km.NN.tif'
output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\PRECIP.csv'

pop.path.and.prefix  <- 'Z:\\GBD2017\\Processing\\Stages\\03_Muster_Population_Figures\\Verified_Outputs\\Ouput_Pop_Unmasked_5K\\ihme_corrected_frankenpop_All_Ages_3_' # the path and all filename characters that precede the year - this will need to parse this according to year, sex, and age
pop.suffix           <- '.tif'  # the filename characters and extension information that follows the year
zone.filename        <- 'Z:\\master_geometries\\Admin_Units\\Global\\MAP\\2018\\MG_5K\\Rasters\\admin2018_1_MG_5K.tif' # this will also be a raster (and will be clipped to the extent of the raster
start.year           <- 2000
end.year             <- 2017
end.replicated.year  <- 2023

n.years <- end.year - start.year + 1

admin <- raster(zone.filename)

template <- raster(template.filename)

for (year in start.year:end.year) {

  if (year == start.year){
    # Establish a template file for subsetting the global rasters (by choosing a realization the whole stack won't need to be subset)
    e <- extent(template)
    r <- res(template)

    # Read and crop in the global rasters
    admin <- crop(admin,e)

    admin.r <- res(admin)
    if ((admin.r[1] != r[1]) || (admin.r[2] != r[2])) {
      admin  <- resample(admin, template, method = 'ngb')
    }

    # Use the limits raster with the admin raster to define the set of endemic countries (only these will be processed on a per country basis - thus reducing runtime)
    endemic.vec         <- unique(admin)
    endemic.vec         <- endemic.vec[!endemic.vec < 1] # Remove the negative and zero values (likely no-data) 
    #omit.list           <- c(1993, 2029, 2042, 2046, 2085, 2086, 2114, 2135, 2160, 2162, 2169, 2174, 2185, 2192, 2204, 2225) # for Africa only - hack addition
    #endemic.vec         <- endemic.vec[!endemic.vec %in% omit.list]
    omit.list           <- c(90, 224, 226, 230, 234, 235, 713, 714, 715, 716, 717, 718, 719, 721, 723, 724, 726, 727, 728, 729, 730, 731, 732, 734, 735, 737, 740, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 1010, 1013, 1018, 1021, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1062, 1149, 1313, 1314, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1439, 1440, 1441, 1442, 1443, 1444, 1445, 1446, 1447, 1448, 1449, 1450, 1451 ,1452, 1582, 1583, 1589, 1592, 1601, 1606, 1614, 1617, 1618, 1620, 1623, 1625, 1629, 1634, 1638, 1639, 1640, 1643, 1645, 1649, 1652, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1950, 1951, 1953, 1959) # for Africa only - hack addition
    endemic.vec         <- endemic.vec[!endemic.vec %in% omit.list]    
    n.endemic.countries <- length(endemic.vec)

    n.output.rows <- n.endemic.countries
    n.output.cols <- 1 + n.years

    # Define the results table - this will contain columns for each realization, but those will be collapsed prior to writing out to only retain the mean, median, LCI, and UCI values
    results.table <- as.data.frame(matrix(NA, nrow = n.output.rows, ncol = n.output.cols))
    
    colnames(results.table)[1] <- 'ID'
  }

  input.filename <- paste(raster.path, year, filename.string, sep='') 
  input <- raster(input.filename)
  input <- crop(input, e)

  pop.filename   <- paste(pop.path.and.prefix, year, pop.suffix, sep='')
  pop <- raster(pop.filename)
  pop   <- crop(pop, e)

  # Check to see for exact pixel alignment and resample if the static files don't perfectly align with the realization
  input.r <- res(input)
  if ((input.r[1] != r[1]) || (input.r[2] != r[2])) {
    input  <- resample(input, input, method = 'ngb')
  }
  pop.r <- res(pop)
  if ((pop.r[1] != r[1]) || (pop.r[2] != r[2])) {
    pop  <- resample(pop, input, method = 'ngb')
  }

  count.raster <- input * pop

  # use the zonal stat shortcut
  full.zonal <- zonal(count.raster, admin, fun='sum')
  pop.zonal  <- zonal(pop, admin, fun='sum')

  rate.zonal <- full.zonal
  rate.zonal[,2] <- full.zonal[,2] / pop.zonal[,2]

  to.keep <- rate.zonal[,1] %in% endemic.vec

  if (year == start.year) { results.table[,1] <- endemic.vec }
  col.pos <- year - start.year + 2 # plus 2 to skip over the ID column
  results.table[,col.pos] <- rate.zonal[to.keep,]

  colnames(results.table)[col.pos] <- toString(year)

  } # end A loop
}

for (year.2 in (end.year + 1) : end.replicated.year) {
  col.pos.2 <- year.2 - start.year + 2 # plus 2 to skip over the ID column
  results.table[,col.pos.2] <- results.table[,col.pos] # borrow the column position from the previous loop
  colnames(results.table)[col.pos.2] <- toString(year.2) 
}

# Write the summary and full realizations tables
write.csv(results.table, output.filename , row.names=F)











### For a single raster
#library(raster)
#template.filename    <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries\\rasters\\PfPR_mean\\pr_2000_rmean_Africa.tif'
#zone.filename        <- 'Z:\\master_geometries\\Admin_Units\\Global\\MAP\\2018\\MG_5K\\Rasters\\admin2018_1_MG_5K.tif' # this will also be a raster (and will be clipped to the extent of the raster
#output.filename      <- 'Z:\\GBD2017\\Processing\\Stages\\09a_Africa_Cubes\\Checkpoint_Outputs\\Summary_20181122\\Africa\\summaries_R_0\\tables\\persistence_results.csv'

#admin <- raster(zone.filename)
#template <- raster(template.filename)

#e <- extent(template)
#r <- res(template)

# Read and crop in the global rasters
#admin <- crop(admin,e)

#admin.r <- res(admin)
#if ((admin.r[1] != r[1]) || (admin.r[2] != r[2])) {
#  admin  <- resample(admin, template, method = 'ngb')
#}

#endemic.vec         <- unique(admin)
#endemic.vec         <- endemic.vec[!endemic.vec < 1] # Remove the negative and zero values (likely no-data) 
#omit.list           <- c(90, 224, 226, 230, 234, 235, 713, 714, 715, 716, 717, 718, 719, 721, 723, 724, 726, 727, 728, 729, 730, 731, 732, 734, 735, 737, 740, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 1010, 1013, 1018, 1021, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1062, 1149, 1313, 1314, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1439, 1440, 1441, 1442, 1443, 1444, 1445, 1446, 1447, 1448, 1449, 1450, 1451 ,1452, 1582, 1583, 1589, 1592, 1601, 1606, 1614, 1617, 1618, 1620, 1623, 1625, 1629, 1634, 1638, 1639, 1640, 1643, 1645, 1649, 1652, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1950, 1951, 1953, 1959) # for Africa only - hack addition
#endemic.vec         <- endemic.vec[!endemic.vec %in% omit.list]    
#n.endemic.countries <- length(endemic.vec)
#n.output.rows <- n.endemic.countries

#results.table <- as.data.frame(matrix(NA, nrow = n.output.rows, ncol = 2))
#colnames(results.table)[1] <- 'ID'

#input.filename <- 'D:\\temp\\novartis\\persistence.tif'
#input <- raster(input.filename)
#input <- crop(input, e)

#pop.filename   <-  'Z:\\GBD2017\\Processing\\Stages\\03_Muster_Population_Figures\\Verified_Outputs\\Ouput_Pop_Unmasked_5K\\ihme_corrected_frankenpop_All_Ages_3_2017.tif'
#pop <- raster(pop.filename)
#pop   <- crop(pop, e)

#count.raster <- input * pop

## use the zonal stat shortcut
#full.zonal <- zonal(count.raster, admin, fun='sum')
#pop.zonal  <- zonal(pop, admin, fun='sum')

#rate.zonal <- full.zonal
#rate.zonal[,2] <- full.zonal[,2] / pop.zonal[,2]

#to.keep <- rate.zonal[,1] %in% endemic.vec

#results.table <- rate.zonal[to.keep,]

## Write the summary and full realizations tables
#write.csv(results.table, output.filename , row.names=F)














