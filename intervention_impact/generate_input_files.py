
import pandas as pd
import pdb
import os
import shapely
import sys
import json

from dtk.tools.demographics.Node import Node, nodeid_from_lat_lon
from simtools.SetupParser import SetupParser
from input_file_generation.DemographicsGenerator import DemographicsGenerator
from input_file_generation.ClimateGenerator import ClimateGenerator
from simtools.Utilities.COMPSUtilities import COMPS_login

sys.path.insert(0, os.path.pardir)
from spatial import make_shapefile, extract_latlongs


SetupParser.init()
COMPS_login("https://comps.idmod.org")

res = 30
pop = 1000
vector_raster_dir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
                       'map_intervention_impact', 'seasonal_classification', 'vectors')

sites = pd.read_csv("site_details.csv")

# # for African sites, find vector species mix
# print("finding vector mix")
# shp_df = make_shapefile(sites.copy(), type='point',
#                             to_crs={'init' :'epsg:4326'}, lat_name='lat', lon_name='lon')
# shapes = [shapely.geometry.mapping(g) for g in shp_df['geometry']]
#
# for species in ['arabiensis', 'funestus', 'gambiae']:
#     sites[species] = extract_latlongs(os.path.join(vector_raster_dir, '{name}.tif'.format(name=species)), shapes)

print("generating input files")
for idx, name in enumerate(sites['name']):

    print("site: " + name)

    # demographics
    print("demographics")
    demog_path = 'demog/demog_{name}.json'.format(name=name)
    site = sites.iloc[idx]
    node_id = nodeid_from_lat_lon(site['lat'], site['lon'], res)
    node = Node(site['lat'], site['lon'], pop, node_id)
    site_demog = DemographicsGenerator([node], res_in_arcsec=res)

    demographics = site_demog.generate_demographics()

    demo_f = open(demog_path, 'w+')
    json.dump(demographics, demo_f, indent=4)
    demo_f.close()

    # climate
    print("climate")
    climate_path = 'climate/{name}'.format(name=name)
    site_climate = ClimateGenerator(demog_path, os.path.join(climate_path, 'climate_wo.json'), climate_path,
                                    climate_project='IDM_{country}'.format(country=site['country']))
    site_climate.generate_climate_files()

    pdb.set_trace()



# DemographicsGenerator.from_grid_file("site_details.csv", "new_demog.json")
