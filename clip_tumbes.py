import os
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import pdb
import seaborn as sns
import shapely.geometry
from sklearn.cluster import KMeans

from spatial import run_all_clipping, plot_all_shapes, make_shapefile, plot_shape, extract_latlongs

rasters_to_clip = {
    "friction":{"input_path": "Z:/mastergrids/Other_Global_Covariates/Accessibility/Weiss/",
              "input_pattern": "friction.*tif$"}
}

main_path = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/peru_northcoast/"
out_path = os.path.join(main_path, "friction_surface/")
shapefile_name = "bbox"
data_path = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/data/Peru/AntonioQuispe_VisitingScholar2018/NorthCoast_Peru/northcoast1.csv"
concave_alpha = 10

print("loading and cleaning point data")
hf_data = pd.read_csv(data_path, encoding='latin1')
hf_data = hf_data.query('year==2011 & week==1') # gives us a unique set of lat-longs

# hf_data.rename(columns={"mid_x": "lng", "mid_y": "lat"}, inplace=True)

for input_name, input_vals in rasters_to_clip.items():

    print("extracting for " + input_name)

    # plot_all_shapes(hf_data, alphas=[1, 5, 10])
    run_all_clipping(input_vals["input_path"], out_path, shapefile_name, hf_data,
                     raster_pattern=input_vals["input_pattern"], overwrite=True,
                     raster_folder=input_name, alpha=concave_alpha, out_name="{name}_all".format(name=input_name),
                     write_shp=True)

    clip_catchments = False
    if clip_catchments:

        for  hf in pd.unique(hf_data['catchment']):
            run_all_clipping(os.path.join(out_path, "rasters", input_name), out_path, shapefile_name,
                             hf_data.query('catchment==@hf'), raster_pattern="{name}_all.*tif$".format(name=input_name),
                             overwrite=True, raster_folder=input_name, alpha=concave_alpha,
                             out_name="{name}_{hf_name}_{shp_type}".format(name=input_name, hf_name=hf, shp_type=shapefile_name),
                             write_shp=False, crop=False)

