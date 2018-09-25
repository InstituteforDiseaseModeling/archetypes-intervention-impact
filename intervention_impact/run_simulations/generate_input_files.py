import numpy as np
import scipy.integrate as sp
import pandas as pd
import pdb
import os
import shapely.geometry
import sys
import re
import json
import shutil

from dtk.tools.demographics.Node import Node, nodeid_from_lat_lon
from simtools.SetupParser import SetupParser
from input_file_generation.DemographicsGenerator import DemographicsGenerator
from input_file_generation.ClimateGenerator import ClimateGenerator
from simtools.Utilities.COMPSUtilities import COMPS_login
from input_file_generation.add_properties_to_demographics import generate_demographics_properties, check_df_valid

sys.path.insert(0, '../..')
from spatial import make_shapefile, extract_latlongs


def update_demog(demographics):

    # impose heterogeneous biting risk:
    demographics["Defaults"]["IndividualAttributes"]["RiskDistributionFlag"] = 3
    demographics["Defaults"]["IndividualAttributes"]["RiskDistribution1"] = 1

    return demographics

def net_usage_overlay(base_demog_path, overlay_path):

    # base IPs and NPs
    IPs = [
        { 'Property' : 'NetUsage',
          'Values' : [ 'HatesNets',
                       'LovesNets'],
          'Initial_Distribution' : [0.1, 0.9],
          'Transitions' : [] }
    ]
    NPs = [
    ]

    generate_demographics_properties(base_demog_path, overlay_path, IPs=IPs, NPs=NPs, as_overlay=True)


def generate_input_files(out_dir, res=30, pop=1000, overwrite=False):

    # setup
    # SetupParser.init()
    COMPS_login("https://comps.idmod.org")

    vector_raster_dir = os.path.join(os.path.expanduser("~"), "Dropbox (IDM)", "Malaria Team Folder", "projects",
                                     "map_intervention_impact", "seasonal_classification", "vectors")

    # specify directories

    if overwrite and os.path.isdir(out_dir):
        shutil.rmtree(out_dir)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    sites = pd.read_csv("site_details.csv")

    print("finding vector mix")
    africa_sites = sites.query("continent=='Africa'")
    shp_df = make_shapefile(africa_sites.copy(), type="point",
                            to_crs={"init": "epsg:4326"}, lat_name="lat", lon_name="lon")
    shapes = [shapely.geometry.mapping(g) for g in shp_df["geometry"]]

    africa_sites = africa_sites[["name"]]
    for species in ["arabiensis", "funestus", "gambiae"]:
        species_prop = extract_latlongs(os.path.join(vector_raster_dir, "{name}.tif".format(name=species)), shapes)
        africa_sites[species] = species_prop

    non_africa_sites = pd.read_csv("vector_props_non_africa.csv")
    all_sites = africa_sites.append(non_africa_sites, sort=False).fillna(0)

    all_sites.to_csv(os.path.join(out_dir, "vector_proportions.csv"), index=False)


    print("generating input files")

    # demographics
    print("demographics")
    demog_path = os.path.join(out_dir, "demographics.json")
    sites['node_id'] = sites.apply(lambda row: nodeid_from_lat_lon(row["lat"], row["lon"], res), axis=1)
    nodes = [Node(this_site["lat"], this_site["lon"], pop,
                  this_site["node_id"], extra_attributes = {"Country": this_site["birth_rate_country"]})
             for ix, this_site in sites.iterrows()]

    site_demog = DemographicsGenerator(nodes, res_in_arcsec=res,
                                       update_demographics=update_demog)
    demographics = site_demog.generate_demographics()

    demo_f = open(demog_path, "w+")
    json.dump(demographics, demo_f, indent=4)
    demo_f.close()

    overlay_path = os.path.join(out_dir, "demographics_net_overlay.json")
    net_usage_overlay(demog_path, overlay_path)