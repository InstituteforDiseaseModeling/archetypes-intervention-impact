import os
import pdb
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import seaborn as sns
import shapely.geometry
from sklearn.cluster import KMeans

from spatial import run_all_clipping, plot_all_shapes, make_shapefile, plot_shape, extract_latlongs

main_path = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/accessibility_example/"
shapefile_name = "bbox"
data_paths = ["seattle.csv", "rainier.csv"]

raster= {"input_path": os.path.join(main_path, "full_friction_surface"),
         "input_pattern": "friction_surface_2015_v1\.0\.tif$"}

data_centers={
    "seattle" : pd.DataFrame({'lat': [47.622650], 'lng': [-122.349195]}),
    "rainier" : pd.DataFrame({'lat': [46.870001], 'lng': [-121.757489]})
}

# for data_name, data in data_centers.items():
for data_name in data_paths:

    print("extracting for " + data_name)
    data = pd.read_csv(os.path.join(main_path, data_name))

    # plot_all_shapes(hh_data, alphas=[1, 5, 10])
    run_all_clipping(raster["input_path"], main_path, shapefile_name, data,
                     raster_pattern=raster["input_pattern"], overwrite=True,
                     out_name="{name}_access".format(name=data_name[:-4]),
                     write_shp=True, alpha=200)
