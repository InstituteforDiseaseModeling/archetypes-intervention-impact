# compare MAP's overall 2015 incidence rate to true incidence and treated cases from an IDM model run

import pandas as pd
import os
import json
import shapely
import rasterio
import pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from simtools.Utilities.COMPSUtilities import COMPS_login
from simtools.Utilities.Experiments import retrieve_simulation
from spatial import make_shapefile, extract_latlongs, mask_raster


print("extracting data from simulation")
COMPS_login("https://comps.idmod.org")

sim_id = "4628a73e-9d39-e811-a2bf-c4346bcb7274"
sim = retrieve_simulation(sim_id)
sim_path = os.path.join(sim.get_path(), "output")
hf = sim.tags["__site__"]

with open(os.path.join(sim_path, "ReportEventCounter.json")) as rec:
    treated_dict = json.loads(rec.read())

with open(os.path.join(sim_path, "ReportMalariaFiltered.json")) as rmf:
    filtered_dict = json.loads(rmf.read())
idm_output = pd.DataFrame({channel: filtered_dict['Channels'][channel]['Data'] for channel in filtered_dict['Channels'].keys()})
idm_output['Received Treatment'] = treated_dict['Channels']['Received_Treatment']['Data']

# these sims were run from 2009-01-01 (day 0) to 2020-12-31 (day 4015). We want 2015, i.e. days 2190 to 2554
idm_output = idm_output.iloc[365*6:365*7]

# plot IDM grids
print("extracting and plotting IDM grid locations and populations")
main_dir = "C:/Users/abertozzivilla/Dropbox (IDM)/Malaria Team Folder/projects/Mozambique/"
raster_dir = os.path.join(main_dir, "incidence_calibration/rasters/")
gridded_inputs = pd.read_csv(os.path.join(main_dir, "gridded_simulation_input/grid_lookup_friction.csv"))
gridded_pop = pd.read_csv(os.path.join(main_dir, 'gridded_simulation_input/grid_population.csv'))
gridded_inputs = gridded_inputs.query('catchment==@hf')
gridded_inputs = pd.merge(gridded_inputs, gridded_pop, how='inner')

fig= plt.figure(figsize=(5,5))
idm_ax = plt.subplot(111)
idm = idm_ax.scatter(gridded_inputs['mid_x'], gridded_inputs['mid_y'],
                      c=gridded_inputs['pop'], vmin=0, marker=',', cmap='viridis')

plt.colorbar(idm, ax=idm_ax)
idm_ax.set_title("IDM Population, Motaze")
idm_ax.xaxis.set_visible(False)
idm_ax.yaxis.set_visible(False)
plt.savefig('idm_pop.png')

# load incidence and population raster for Motaze
print("loading MAP rasters")
inc_dir = os.path.join(raster_dir, "incidence/incidence_all.year.2015.tif")
pop_dir = os.path.join(raster_dir, "pop_5k/pop_5k_all.year.2015.tif")
map_inc_raster = rasterio.open(inc_dir)
map_pop_raster = rasterio.open(pop_dir)

# mask to idm grid cells only
print("masking map rasters")
sim_shp = make_shapefile(gridded_inputs.copy(), type='point',
                            to_crs={'init' :'epsg:4326'},
                            lat_name='mid_y', lon_name='mid_x')
shapes = [shapely.geometry.mapping(g) for g in sim_shp['geometry']]

map_inc, map_inc_meta = mask_raster(map_inc_raster, shapes)
map_pop, clip_pop_meta = mask_raster(map_pop_raster, shapes)

map_inc = map_inc.data[0]
map_pop = map_pop.data[0]


# all cases
print("plotting masked rasters")
map_count_cases = pd.DataFrame({'inc': map_inc[map_inc>0],
                               'pop': map_pop[map_pop>0]})
map_count_cases['cases'] = map_count_cases['inc'] * map_count_cases['pop']

fig= plt.figure(figsize=(10,5))
inc_ax = plt.subplot(121)
pop_ax = plt.subplot(122)

masked_array = np.ma.array(map_inc, mask=map_inc<0)
cmap = mpl.cm.viridis
cmap.set_bad('white',1.)
inc_plot = inc_ax.imshow(masked_array, vmin=0, cmap=cmap)
plt.colorbar(inc_plot, ax=inc_ax)
inc_ax.set_title("MAP Incidence, Motaze")
inc_ax.xaxis.set_visible(False)
inc_ax.yaxis.set_visible(False)

masked_array = np.ma.array(map_pop, mask=map_pop==0)
pop_plot = pop_ax.imshow(masked_array, vmin=0, cmap=cmap)
plt.colorbar(pop_plot, ax=pop_ax)
pop_ax.set_title("MAP Population, Motaze")
pop_ax.xaxis.set_visible(False)
pop_ax.yaxis.set_visible(False)

plt.savefig('map_inc_pop_point_new.png')

# load true cases
print("loading observed cases")
catch_inc = pd.read_csv(os.path.join(main_dir, "gridded_simulation_input/catchment_incidence.csv"))
catch_inc = catch_inc.query('catchment==@hf')
catch_inc['year'] = catch_inc['fulldate'].apply(lambda x: int(x[0:4]))

# find cases and population-normalized cases for each

summary = {'data': {'cases': catch_inc.query('year==2015')['cases'].sum(),
                    'pop': gridded_inputs['pop'].sum()},
           'map': {'cases': map_count_cases.query('pop<2000')['cases'].sum(),
                   'pop': map_count_cases.query('pop<2000')['pop'].sum()},
           'idm_new_infections': {'cases': idm_output['New Infections'].sum(),
                                  'pop': idm_output['Statistical Population'].mean()},
           'idm_clinical_cases': {'cases': idm_output['New Clinical Cases'].sum(),
                                  'pop': idm_output['Statistical Population'].mean()},
           'idm_received_treatment': {'cases': idm_output['Received Treatment'].sum(),
                                      'pop': idm_output['Statistical Population'].mean()}
           }

summary = pd.DataFrame.from_dict(summary, orient='index')
summary['incidence_per_1000'] = summary['cases']/summary['pop'] *1000
summary['pop_normalized_cases'] = summary['incidence_per_1000']*summary.loc['data']['pop']/1000

print("Results:")
print(summary)
summary.to_csv('inc_map_summary.csv')

