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

def wavg(group, avg_name, weight_name):
    """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
    In rare instance, we may not have weights, so just return the mean. Customize this if your business case
    should return otherwise.
    """
    d = group[avg_name]
    w = group[weight_name]
    try:
        return (d * w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()

main_dir = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/Mozambique/"
raster_dir = os.path.join(main_dir, "incidence_calibration/rasters/")
out_dir = os.path.join(main_dir, "incidence_calibration/transmission_population_comparison/")
gridded_inputs = pd.read_csv(os.path.join(main_dir, "gridded_simulation_input/grid_lookup_friction.csv"))
gridded_inputs = gridded_inputs.query('catchment!="Caputine"') #todo: is this my bug?

# load gold standard data
pop_data = pd.read_csv(os.path.join(main_dir, 'gridded_simulation_input/grid_population.csv'))
pop_data = pd.merge(pop_data, gridded_inputs, how='left')
grid_counts = pop_data.groupby('catchment').count()[['grid_cell']]
sim_shp = make_shapefile(gridded_inputs.copy(), type='point',
                         to_crs={'init': 'epsg:4326'},
                         lat_name='mid_y', lon_name='mid_x')
shapes = [shapely.geometry.mapping(g) for g in sim_shp['geometry']]

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
prev_data = pd.merge(prev_data.reset_index(), pop_data)
prev_data['type'] = "data"
catch_data_vals = prev_data.groupby(['type', 'catchment']).sum()[['sum', 'count']].reset_index()
catch_data_vals['prev'] = catch_data_vals['sum']/catch_data_vals['count']
catch_data_vals.drop(['sum', 'count'], axis=1, inplace=True)
catch_data_pop = pop_data.groupby('catchment').sum()['pop'].reset_index()
catch_data_pop.rename(columns={'pop': 'census_pop'}, inplace=True)


print("extracting and summarizing custom Mozambique surface")
ewan_vals = pd.read_csv(os.path.join(raster_dir, "ewan_prev/grid_vals.csv"))
ewan_vals = pd.pivot_table(ewan_vals, index=['grid_cell', 'catchment'], values='ewan_prev', columns=['type']).reset_index()
ewan_vals = ewan_vals.query('mean>0')
ewan_vals['ewan_pop'] = ewan_vals['popdens']/4.6**2
ewan_vals['lower'] = ewan_vals['mean'] - ewan_vals['sd']*1.96
ewan_vals['upper'] = ewan_vals['mean'] + ewan_vals['sd']*1.96
ewan_vals = ewan_vals.melt(id_vars=['catchment', 'grid_cell', 'ewan_pop'], value_vars=['mean', 'lower', 'upper'], value_name="prev")
catch_ewan_vals = ewan_vals.groupby(['catchment', 'type']).apply(wavg, "prev", "ewan_pop").reset_index().rename(columns={0:"prev"})
catch_ewan_pop = ewan_vals.query('type=="mean"').groupby('catchment').sum()['ewan_pop'].reset_index()
catch_ewan_vals['type'] = catch_ewan_vals['type'].apply(lambda x: "ewan_"+x)

print("extracting and summarizing MAP surface")
map_vals = pd.read_csv(os.path.join(raster_dir, "prevalence/grid_vals.csv"))
map_pop =  pd.read_csv(os.path.join(raster_dir, "pop_5k/grid_vals.csv"))

map_vals = pd.merge(map_vals, map_pop).query('prevalence>0')
map_vals['afri_pop'] = map_vals['pop_5k']/25
catch_map_vals = map_vals.groupby(['sample', 'catchment']).apply(wavg, "prevalence", "afri_pop").reset_index().rename(columns={0:"prev"})
catch_map_pop = map_vals.query("sample==1").groupby('catchment').sum()['afri_pop'].reset_index()
catch_map_vals['type'] = catch_map_vals['sample'].apply(lambda x: "map_sample_" + str(x))
catch_map_vals.drop('sample', axis=1, inplace=True)

print("cleaning FB pop")
fb_pop = pd.read_csv(os.path.join(out_dir, "fb_pop.csv"))
fb_pop =  fb_pop.dropna().groupby('catchment').sum()['fb_pop'].reset_index()

print("joining prevalence datasets")
all_prev = pd.concat([catch_data_vals, catch_ewan_vals, catch_map_vals])
all_prev = all_prev[["type", "catchment", "prev"]]

print("joining population datasets")
all_pop = pd.merge(pd.merge(pd.merge(catch_data_pop, catch_ewan_pop), catch_map_pop), fb_pop)

full_dataset = pd.merge(all_prev, all_pop)
full_dataset.to_csv(os.path.join(out_dir, "full_catchment_compare.csv"), index=False)