import os
import pandas as pd
import re
import shapely
import rasterio
import pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime
from dateutil.relativedelta import relativedelta
import seaborn as sns

from spatial import make_shapefile, extract_latlongs, mask_raster


print("extracting and plotting IDM grid locations and populations")
main_dir = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/Mozambique/"
raster_dir = os.path.join(main_dir, "incidence_calibration/rasters/")
gridded_inputs = pd.read_csv(os.path.join(main_dir, "gridded_simulation_input/grid_lookup_friction.csv"))
gridded_inputs = gridded_inputs.query('catchment!="Caputine"') #todo: is this my bug?

raster_dirs = {"synoptic_temp": "month"}

# shared datasets
pop_data = pd.read_csv(os.path.join(main_dir, 'gridded_simulation_input/grid_population.csv'))
pop_data = pd.merge(pop_data, gridded_inputs, how='left')
grid_counts = pop_data.groupby('catchment').count()[['grid_cell']]
sim_shp = make_shapefile(gridded_inputs.copy(), type='point',
                         to_crs={'init': 'epsg:4326'},
                         lat_name='mid_y', lon_name='mid_x')
shapes = [shapely.geometry.mapping(g) for g in sim_shp['geometry']]

for dir, unit in raster_dirs.items():
    print("extracting for " + dir)

    files = os.listdir(os.path.join(raster_dir, dir))
    pattern = re.compile(".*\.tif$")
    tifs = [x for x in files if pattern.match(x)]
    tifs.sort()

    cov_list = []

    for idx, tif in enumerate(tifs):
        print(tif)
        this_df = gridded_inputs.copy()
        this_df["{name}".format(name=dir)] = extract_latlongs(os.path.join(raster_dir, dir, tif), shapes)
        this_df[unit] = idx + 1

        cov_list.append(this_df)

    cov_df = pd.concat(cov_list)
    cov_df.to_csv(os.path.join(raster_dir, dir, "extracted_vals.csv"), index=False)


