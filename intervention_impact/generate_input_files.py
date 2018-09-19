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

sys.path.insert(0, os.path.pardir)
from spatial import make_shapefile, extract_latlongs


def update_demog(demographics):

    # impose heterogeneous biting risk:
    demographics["Defaults"]["IndividualAttributes"]["RiskDistributionFlag"] = 3
    demographics["Defaults"]["IndividualAttributes"]["RiskDistribution1"] = 1

    return demographics

def net_usage_overlay(site_name, hates_net_prop):

    base_demog_path = os.path.join("sites", site_name, "demographics_{name}.json".format(name=site_name))
    overlay_path = os.path.join("sites", site_name, "demographics_{name}_hatenets_{prop}.json".format(name=site_name,
                                                                                                        prop=hates_net_prop))
    with open(base_demog_path) as f:
        demo = json.loads(f.read())

    nodeid = demo["Nodes"][0]["NodeID"]

    prop_by_node_dict = { 'node' : [nodeid, nodeid],
                          'Property' : ['NetUsage', 'NetUsage'],
                          'Property_Type' : ['IP', 'IP'],
                          'Property_Value' : [ 'HatesNets',
                                               'LovesNets'],
                          'Initial_Distribution' : [hates_net_prop, 1-hates_net_prop]}

    # base IPs and NPs
    IPs = [
        { 'Property' : 'NetUsage',
          'Values' : [ 'HatesNets',
                       'LovesNets'],
          'Initial_Distribution' : [0, 1],
          'Transitions' : [] }
    ]
    NPs = [
    ]

    df = pd.DataFrame(prop_by_node_dict)
    if (not df.empty) and (not check_df_valid(df, IPs, NPs)):
        print('properties by node df is invalid')
        exit()
    generate_demographics_properties(base_demog_path, overlay_path, IPs=IPs, NPs=NPs, df=df, as_overlay=True)



def generate_input_files(site_name, res=30, pop=1000, overwrite=False):

    # setup
    # SetupParser.init()
    COMPS_login("https://comps.idmod.org")

    vector_raster_dir = os.path.join(os.path.expanduser("~"), "Dropbox (IDM)", "Malaria Team Folder", "projects",
                                     "map_intervention_impact", "seasonal_classification", "vectors")

    # specify directories
    out_dir = os.path.join("sites", site_name)

    if overwrite and os.path.isdir(out_dir):
        shutil.rmtree(out_dir)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    sites = pd.read_csv("site_details.csv")
    this_site = sites.query("name==@site_name")

    print("finding vector mix")
    if this_site["continent"].iloc[0] == "Africa":
        shp_df = make_shapefile(this_site.copy(), type="point",
                                to_crs={"init": "epsg:4326"}, lat_name="lat", lon_name="lon")
        shapes = [shapely.geometry.mapping(g) for g in shp_df["geometry"]]

        site_vectors = pd.DataFrame(columns=["species", "proportion"])
        for species in ["arabiensis", "funestus", "gambiae"]:
            species_prop = extract_latlongs(os.path.join(vector_raster_dir, "{name}.tif".format(name=species)), shapes)[0]
            if species_prop > 0:
                site_vectors = site_vectors.append({"species": species,
                                                    "proportion": species_prop},
                                                     ignore_index=True)
    else:
        site_vectors = pd.read_csv("vector_props_non_africa.csv")
        site_vectors = site_vectors[["species", site_name]].rename(columns={site_name: "proportion"})
        site_vectors = site_vectors.query("proportion>0")

    site_vectors.to_csv(os.path.join(out_dir, "vector_proportions.csv"), index=False)

    print("generating input files")

    # demographics
    print("demographics")
    demog_path = os.path.join(out_dir, "demographics_{name}.json".format(name=site_name))
    node_id = nodeid_from_lat_lon(this_site["lat"], this_site["lon"], res)
    node = Node(this_site["lat"], this_site["lon"], pop, node_id)
    site_demog = DemographicsGenerator([node], res_in_arcsec=res, country=this_site["birth_rate_country"].iloc[0],
                                       update_demographics=update_demog)
    demographics = site_demog.generate_demographics()

    demo_f = open(demog_path, "w+")
    json.dump(demographics, demo_f, indent=4)
    demo_f.close()

    # climate
    print("climate")
    orig_country_name = this_site["country"].iloc[0]
    site_climate = ClimateGenerator(demog_path, os.path.join(out_dir, "climate_wo.json"), out_dir,
                                    climate_project="IDM-{country}".format(country=orig_country_name),
                                    start_year=2016)
    site_climate.generate_climate_files()

    # rename climate files
    print("renaming climate files")

    # tedious: Name of output files not consistent with input name, need to manually change
    out_country_names = {"Burkina_Faso": "Burkina Faso",
                         "Democratic_Republic_of_the_Congo": "DemocraticRepublicOfTheCongo",}

    country = out_country_names[orig_country_name] if orig_country_name in out_country_names.keys() else orig_country_name

    res_unit = "arcsec" if res == 30 else "arcmin"
    for fname in os.listdir(out_dir):
        match = "{country}_{res}{res_name}_(.*)".format(country=country, res=res, res_name=res_unit)
        new_name = re.sub(match, r"\1", fname)
        os.replace(os.path.join(out_dir, fname),
                   os.path.join(out_dir, new_name))

        # delete comps-generated demog file
        if "arcmin_demographics" in fname:
            os.remove(os.path.join(out_dir, fname))
