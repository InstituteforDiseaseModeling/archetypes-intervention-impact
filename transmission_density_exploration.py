import pandas as pd
import os
import json
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

print("loading MAP rasters")
inc_dir = os.path.join(raster_dir, "incidence/incidence_all.year.2015.tif")
pop_dir = os.path.join(raster_dir, "pop_5k/pop_5k_all.year.2015.tif")
inc_raster_dir = os.path.join(raster_dir, "baseline_incidence")
prev_raster_dir = os.path.join(raster_dir, "baseline_prevalence")

print("incidence analysis")
# load catchment incidence
print("extracting data")
inc_data = pd.read_csv(os.path.join(main_dir, 'gridded_simulation_input/catchment_incidence_error_correction.csv'))
pop_data = pd.read_csv(os.path.join(main_dir, 'gridded_simulation_input/grid_population.csv'))
pop_data = pd.merge(pop_data, gridded_inputs, how='left')
grid_counts = pop_data.groupby('catchment').count()[['grid_cell']]
catch_pop_data = pop_data.groupby('catchment').sum()[['pop']]
catch_pop_data = pd.merge(catch_pop_data, grid_counts, left_index=True, right_index=True)
inc_data['year'] = inc_data['fulldate'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d').year)
inc_data.drop('fulldate', axis=1, inplace=True)
inc_data.loc[inc_data.catchment=='Panjane-Caputine', 'catchment'] = "Caputine-Panjane"
inc_data = inc_data.groupby(['year', 'catchment']).sum().loc[2015]
inc_data = pd.merge(inc_data, catch_pop_data, left_index=True, right_index=True)
inc_data['inc_rate'] = inc_data['cases']/inc_data['pop']*1000
inc_data['avg_pop'] = inc_data['pop']/inc_data['grid_cell']

sim_shp = make_shapefile(gridded_inputs.copy(), type='point',
                            to_crs={'init' :'epsg:4326'},
                            lat_name='mid_y', lon_name='mid_x')
shapes = [shapely.geometry.mapping(g) for g in sim_shp['geometry']]
gridded_inputs["map_pop"] = [x/25 for x in extract_latlongs(pop_dir, shapes)]

print("extracting MAP rasters")
inc_files = os.listdir(inc_raster_dir)
pattern = re.compile("baseline_incidence_all.*\.tif$")
inc_tifs = [x for x in inc_files if pattern.match(x)]
inc_tifs.sort()

map_inc_list = []

for idx, tif in enumerate(inc_tifs):
    print(tif)
    this_df = gridded_inputs.copy()
    this_df["map_inc_2015"] = extract_latlongs(os.path.join(inc_raster_dir, tif), shapes)
    this_df['map_cases'] = this_df['map_inc_2015'] * this_df['map_pop']
    if "year" in tif:
        this_df['type'] = "MAP: Mean Across Samples"
    else:
        this_df['type'] = "MAP: Sample {val}".format(val=idx+1)

    map_inc_list.append(this_df)

map_inc_df = pd.concat(map_inc_list)
map_inc_df = map_inc_df.groupby(['type', 'catchment']).sum()[['map_pop', 'map_cases']]
map_inc_df['map_inc_2015'] = map_inc_df['map_cases']/map_inc_df['map_pop']*1000

map_inc_df = pd.merge(map_inc_df, grid_counts, left_index=True, right_index=True)
map_inc_df['map_pop_avg'] = map_inc_df['map_pop']/map_inc_df['grid_cell']

print("plotting for incidence")
#scale limits overall
xmin = 0
xmax = max(map_inc_df['map_inc_2015'].max(), inc_data['inc_rate'].max())
ymin = 0
ymax = max(map_inc_df['map_pop_avg'].max(), inc_data['avg_pop'].max()) + 10

# fig = plt.figure(figsize=(12,6))

fig, axes = plt.subplots(2, int(np.ceil((len(inc_tifs)+1) / 2)), sharex=True, sharey=True, figsize=(15, 10))
axes_list = [item for sublist in axes for item in sublist]

data_ax=axes_list.pop(0)
data_plot = data_ax.scatter(inc_data['inc_rate'], inc_data['avg_pop'],
                      c=inc_data['inc_rate'], vmin=xmin, vmax=xmax, cmap='viridis_r')

data_ax.set_title("Data")
data_ax.set_xlabel("Data: All-Age Incidence")
data_ax.set_ylabel("Average Population, 1k by 1k pixel")

for map_type in map_inc_df.index.levels[0]:

    subset = map_inc_df.loc[map_type]
    map_ax = axes_list.pop(0)
    map = map_ax.scatter(subset['map_inc_2015'], subset['map_pop_avg'],
                          c=subset['map_inc_2015'], vmin=xmin, vmax=xmax, cmap='viridis_r')
    map_ax.set_title(map_type)
    map_ax.set_xlabel("MAP Incidence 2015")
    map_ax.set_ylabel("Afripop, Averaged 1k by 1k Pixel")

final_ax = axes_list.pop(0)
plt.colorbar(data_plot, ax=final_ax)
fig.tight_layout()
plt.savefig("map_inc_analysis_samples.png")
#
#
# # prevalence
print("prevalence analysis")
print("extracting and subsetting data")
# load microdata, extract ages 2-10 for 2015 survey
hh_grid_map = pd.read_csv(os.path.join(main_dir, 'gridded_simulation_input/grid_household_lookup.csv'))
hh_grid_map.rename(columns={'ID': 'ID_hh'}, inplace=True)
prev_data = pd.read_csv(os.path.join(main_dir, '../../data/Mozambique/Magude/Members/census_mda_members.csv'),
                        low_memory=False)
prev_data = prev_data[['tested_rdt_mda1', 'ID', 'ID_hh',
                       'DOB', 'visit_date_mda1', 'rdt_positive_mda1']]
prev_data = prev_data.query('tested_rdt_mda1==1')
prev_data.dropna(inplace=True)
prev_data['DOB'] = prev_data['DOB'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
prev_data['visit_date_mda1'] = prev_data['visit_date_mda1'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
prev_data['age'] = prev_data.apply(lambda row: relativedelta(row['visit_date_mda1'], row['DOB']).years, axis=1)
prev_data = prev_data.query('age>=2 & age<10')
prev_data = pd.merge(prev_data, hh_grid_map, how='left')
prev_data = prev_data.groupby('grid_cell').agg({'rdt_positive_mda1': ['sum', 'count']})['rdt_positive_mda1']
prev_data['prev'] = prev_data['sum']/prev_data['count']
prev_data = pd.merge(prev_data.reset_index(), pop_data)

print("extracting MAP rasters")
prev_files = os.listdir(prev_raster_dir)
pattern = re.compile("baseline_prevalence_all.*\.tif$")
prev_tifs = [x for x in prev_files if pattern.match(x)]
prev_tifs.sort()

map_prev_list = []

for idx, tif in enumerate(prev_tifs):
    print(tif)
    this_df = gridded_inputs.copy()
    this_df["map_prev_2015"] = extract_latlongs(os.path.join(prev_raster_dir, tif), shapes)
    if "year" in tif:
        this_df['type'] = "MAP: Mean Across Samples"
    else:
        this_df['type'] = "MAP: Sample {val}".format(val=idx+1)

    map_prev_list.append(this_df)

map_prev_df = pd.concat(map_prev_list)
map_prev_df = pd.merge(map_prev_df, pop_data)

print("plotting prevalence")

fig, axes = plt.subplots(len(inc_tifs)+1, 2, figsize=(20, 10))
axes_list = [item for sublist in axes for item in sublist]

data_ax=axes_list.pop(0)
data_plot = data_ax.scatter(prev_data['prev'], prev_data['pop'],
                      c=prev_data['prev'], vmin=0, vmax=1, cmap='viridis_r')

data_ax.set_title("Data")
data_ax.set_xlabel("Data: Prevalence, Ages 2-10")
data_ax.set_ylabel("Population")

spatial_ax = axes_list.pop(0)
data_spatial = spatial_ax.scatter(prev_data['mid_x'], prev_data['mid_y'],
                      c=prev_data['prev'], vmin=0, vmax=1, marker=',', cmap='viridis_r')

plt.colorbar(data_spatial, ax=spatial_ax)
spatial_ax.set_title("Data Prevalence, 2015")
spatial_ax.xaxis.set_visible(False)
spatial_ax.yaxis.set_visible(False)

for map_type in map_prev_df['type'].unique():

    subset = map_prev_df.query('type==@map_type')
    map_ax = axes_list.pop(0)

    # map = sns.pairplot(x_vars=['map_prev_2015'], y_vars=['pop'], data=subset, hue="catchment")

    map = map_ax.scatter(subset['map_prev_2015'], subset['pop'],
                          c=subset['map_prev_2015'], vmin=0, vmax=1, cmap='viridis_r')
    map_ax.set_xlim(0, 1)
    map_ax.set_ylim(0, prev_data['pop'].max())
    map_ax.set_title(map_type)
    map_ax.set_xlabel("MAP Prevalence, Ages 2-10")
    map_ax.set_ylabel("Afripop")

    spatial_ax = axes_list.pop(0)
    map_spatial = spatial_ax.scatter(subset['mid_x'], subset['mid_y'],
                                      c=subset['map_prev_2015'], vmin=0, vmax=1, marker=',', cmap='viridis_r')

    plt.colorbar(map_spatial, ax=spatial_ax)
    spatial_ax.set_title(map_type)
    spatial_ax.xaxis.set_visible(False)
    spatial_ax.yaxis.set_visible(False)

fig.tight_layout()
plt.savefig("map_prev_analysis_samples.png")

prev_data.to_csv('prev_data_2_10.csv', index=False)
map_prev_df.to_csv('map_prev.csv', index=False)

