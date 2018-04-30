
import os
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio.plot
import seaborn as sns
import shapely

from covariate_prep.spatial import mask_raster

main_path = "Z:/amelia/larval_habitat/karen/"
grid_path = os.path.join(main_path, "rasters/land_cover")
shp_path = os.path.join(main_path, "shp/point.shp")
data_path = "Q:/Malaria/Myanmar/Karen/metf_data/"

print("loading shapefile")
karen_points = gpd.read_file(shp_path)
shapes = [shapely.geometry.mapping(g) for g in karen_points['geometry']]

# for tsi
#
tsi_path = os.path.join(main_path, "rasters/tsi")
tifs = os.listdir(tsi_path)
tifs.sort()

print("extracting points")
for idx, fname in enumerate(tifs):

    in_raster = rasterio.open(os.path.join(tsi_path, fname))

    month = int(re.match('.*\.([0-9]{2})\.', fname).groups()[0])
    print(month)

    # get raster values at shapefile points
    latlongs = [point['coordinates'] for point in shapes]
    raster_values = in_raster.sample(latlongs)
    values_list = [value[0] for value in raster_values]
    karen_points['tsi_{month}'.format(month=month)] = values_list

karen_points['latitude'] = karen_points['geometry'].apply(lambda geom: geom.x)
karen_points.sort_values(by=['latitude'], inplace=True)
to_plot = karen_points.filter(regex="tsi", axis=1)

fig, ax = plt.subplots(1,1, figsize=(7,7))
heat = sns.heatmap(to_plot, cmap="viridis", ax=ax)
ax.set_title("Synoptic TSI, Karen Villages North to South")
plt.savefig(os.path.join(tsi_path, "tsi_heat.jpg"))

print("loading case data")
inc = pd.read_csv(os.path.join(data_path, "incidence_data.csv"))
prev = pd.read_csv(os.path.join(data_path, "prevalence_data.csv"), encoding='mac_roman')

print("cleaning case data")
inc["year"] = inc["d_start"].apply(lambda x: int(x[-4:]))
summed = inc.groupby(['year', 'mp_code']).sum()
summed = summed.filter(regex="rdt_pf", axis=1)
summed["rdt_all"] = summed.sum(axis=1)

# find population
def find_pop(x):
    if not np.isnan(x[0]):
        return x[0]
    else:
        return x[1]*5
prev['pop'] = prev[['mean_pop_month', 'number_hou']].apply(find_pop, axis=1)
prev = prev[(pd.notnull(prev['latitude'])) & (prev['area_zone']<200)]

inc_rates = pd.merge(summed[['rdt_all']].reset_index(), prev[['mp_code', 'pop']], how='inner')
inc_rates['inc_rate'] = inc_rates['rdt_all'] / inc_rates['pop'] * 1000
inc_rates.rename(columns={'mp_code':'id'}, inplace=True)

print("reading raster")
in_raster = rasterio.open(os.path.join(grid_path, "A2013001_MCD12Q1.tif"))
lc_all = in_raster.read(1)
lc_vals = np.unique(lc_all)
lc_vals = lc_vals[np.where(lc_vals != in_raster.meta['nodata'])].tolist()
lc_vals = [x for x in lc_vals if x not in (6, 10, 11)] # low_information land cover types

fig, ax = plt.subplots(1,1, figsize=(7, 9))
rasterio.plot.show(in_raster, ax=ax)
ax.set_title("Land Cover, Karen 2013")
plt.savefig(os.path.join(grid_path, 'lc_2013.jpg'))

def land_type_prop(in_raster, shape, lc_type):
    shapes = [shapely.geometry.mapping(shape)]

    mask, meta = mask_raster(in_raster, shapes)
    data = mask.data[0].flatten()
    data = data[np.where(data != meta['nodata'])]

    prop = len(np.where(data == lc_type)[0]) / len(data)

    return prop

prop_df_list =  []

for radius in [0.01, 0.02, 0.03, 0.04]:

    print("finding buffers for radius ", str(radius))

    print("loading shapefile")
    karen_points = gpd.read_file(shp_path)
    shapes = [shapely.geometry.mapping(g) for g in karen_points['geometry']]

    print("finding buffers")
    karen_points["buffer"] = karen_points["geometry"].apply(lambda x: x.buffer(radius))

    # karen_buffer = karen_points[['id', 'buffer']]
    # karen_buffer.rename(columns={'buffer': 'geometry'}, inplace=True)
    # karen_buffer = gpd.GeoDataFrame(karen_buffer)
    #
    # fig, ax = plt.subplots(1,1)
    # karen_buffer.plot(ax=ax, alpha=0.5)
    # karen_points.plot(ax=ax, color='black')

    print("finding lc props")
    for lc_type in lc_vals:
        karen_points["lc_{cov}_prop".format(cov=lc_type)] = karen_points["buffer"].apply(lambda x: land_type_prop(in_raster, x, lc_type))

    prop_df = pd.merge(inc_rates[['year', 'id', 'inc_rate']], karen_points.filter(regex="lc|id", axis=1), how='inner')
    prop_df['radius'] = radius
    prop_df = pd.melt(prop_df, id_vars=['radius', 'year', 'id', 'inc_rate'], value_name="prop")
    prop_df['variable'] = prop_df['variable'].apply(lambda x: x[:-5])

    prop_df_list.append(prop_df)

# combine case data and land cover data for regression
for_regression = pd.concat(prop_df_list)

# take year with most villages
time_subset = for_regression.query('year==2016')

g = sns.FacetGrid(time_subset, row='radius', col="variable", margin_titles=True)
g.map(plt.scatter,'prop', 'inc_rate', edgecolor='w')
g.savefig(os.path.join(grid_path, 'lc_props_vs_inc.jpg'))

# remove land cover types that are entirely zero
# time_subset = time_subset.query('variable in ["lc_2", "lc_5", "lc_8", "lc_14"]')
time_subset['ln_inc_rate'] = time_subset['inc_rate'].apply(lambda x: np.log(x))

g = sns.FacetGrid(time_subset, row='radius', col="variable", margin_titles=True)
g.map(plt.scatter,'prop', 'ln_inc_rate', edgecolor='w')
g.savefig(os.path.join(grid_path, 'lc_props_vs_log_inc.jpg'))


# print("reshaping")
# karen_melted = pd.melt(karen_points, id_vars=['id','geometry'], var_name="month", value_name="tsi")
# karen_melted["month"] = karen_melted["month"].apply(lambda str: int(re.sub("tsi_", "", str)))
#
