###############################################################################################################
## 01_generate_demographics.py
## Amelia Bertozzi-Villa
## October 2019
##
## EMOD needs climate files and demographics files. Based on the instructions specified in the input_params.json,
## this script creates a new demographics file with the appropriate vector and population parameters.
##
## Requirements: input directory containing an "input_params.json" and a "site_details.csv" file.
##                see "example_input_dir" in this repo.
##
## Outputs: directories named "demog" and "vector" containing the specified input files.
##############################################################################################################

## Importing and setup ---------------------------------------------------------------------------------------
import pandas as pd
import os
import shapely.geometry
import sys
import re
import json

from dtk.tools.demographics.Node import Node
from input_file_generation.DemographicsGenerator import DemographicsGenerator
from input_file_generation.add_properties_to_demographics import generate_demographics_properties

sys.path.insert(0, '../..')  # add the path to spatial.py
from spatial import make_shapefile, extract_latlongs

desired_width = 320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)

## Functions ---------------------------------------------------------------------------------------

def find_vector_props_africa(africa_sites, vector_raster_dir):
    # Find relative vector abundance of three african vectors
    shp_df = make_shapefile(africa_sites.copy(), type="point",
                            to_crs={"init": "epsg:4326"}, lat_name="lat", lon_name="lon")
    shapes = [shapely.geometry.mapping(g) for g in shp_df["geometry"]]

    africa_props = africa_sites[["cluster"]].copy()
    for species in ["arabiensis", "funestus", "gambiae"]:
        pattern = re.compile(".*_{species}.*\.tif$".format(species=species))
        species_rasters = [x for x in os.listdir(vector_raster_dir) if pattern.match(x)]
        if len(species_rasters)>1:
            print("WARNING: MORE THAN ONE POSSIBLE SPECIES INPUT FILE")

        species_prop = extract_latlongs(os.path.join(vector_raster_dir, species_rasters[0]), shapes)
        africa_props[species] = species_prop

    return africa_props


def find_vector_props_non_africa(non_africa_sites):
    col_names = ["cluster", "darlingi", "minimus", "maculatus"]
    non_africa_vectors = pd.DataFrame(columns=col_names)
    for idx in list(range(0, len(non_africa_sites))):
        site = non_africa_sites.iloc[idx]
        if site["continent"] == "Asia":
            non_africa_vectors = non_africa_vectors.append(pd.DataFrame([[ site["cluster"], 0, 0.6, 0.4]],
                                                                        columns=col_names))
        elif site["continent"] == "Americas":
            non_africa_vectors = non_africa_vectors.append(pd.DataFrame([[site["cluster"], 1, 0, 0]],
                                                                        columns=col_names))
        else:
            print("Continent {cont} not found!".format(cont=site["continent"]))

    return non_africa_vectors


def update_demog(demographics, vectors, tot_vector_count=20000):

    # impose heterogeneous biting risk:
    demographics["Defaults"]["IndividualAttributes"]["RiskDistributionFlag"] = 3
    demographics["Defaults"]["IndividualAttributes"]["RiskDistribution1"] = 1

    # add node-specific vectors
    vectors = pd.melt(vectors, id_vars=["cluster", "continent", "lat", "lon"],
                      var_name="species", value_name="proportion")
    vectors["count"] = vectors["proportion"] * tot_vector_count

    for node in demographics["Nodes"]:
        node_id = node["NodeID"] # used in the pd.query command below
        node["NodeAttributes"]["InitialVectorsPerSpecies"] = {row["species"]: round(row["count"])
                                                                  for idx, row in
                                                                  vectors.query(
                                                                      "cluster==@node_id").iterrows()}
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


## Main code setup ---------------------------------------------------------------------------------------

main_dir = os.path.join(os.path.expanduser("~"),
                        "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                        "20191009_megatrends_era5_new_archetypes")

# uncomment this to run the example
# main_dir = "example_input_dir"

with open(os.path.join(main_dir, "input_params.json")) as f:
    instructions = json.loads(f.read())

main_input_dir = instructions["root_dir"]

sites = pd.read_csv(os.path.join(main_dir, instructions["site_fname"]))

## Demographics ---------------------------------------------------------------------------------------
## Generate main file with site-specific species and equilibrium population dynamics. also add net overlays.

demog_out_dir = os.path.join(main_dir, "demog")
demog_out_path = os.path.join(demog_out_dir, "demographics.json")

if not os.path.isdir(demog_out_dir):
    os.mkdir(demog_out_dir)

if os.path.exists(demog_out_path) and instructions["overwrite_input_files"] == "False":
    print("Demographics files already generated")
else:
    print("Generating Demographics")

    # Step 1: relative vector abundance. Define manually outside of Africa, and using rasters within Africa.
    print("Finding vector mix")
    vector_props = find_vector_props_africa(sites.query("continent=='Africa'"),
                                    os.path.join(os.path.expanduser("~"), instructions["vector_raster_dir"]))

    non_africa = sites.query("continent!='Africa'")
    if len(non_africa)>0:
        non_africa_vectors = find_vector_props_non_africa(non_africa)
        vector_props = vector_props.append(non_africa_vectors, sort=False).fillna(0)

    vector_props = pd.merge(sites[["cluster", "continent", "lat", "lon"]], vector_props)
    vector_out_dir = os.path.join(main_dir, "vector")
    if not os.path.isdir(vector_out_dir):
        os.mkdir(vector_out_dir)
    vector_props.to_csv(os.path.join(vector_out_dir, "vector_proportions.csv"), index=False)

    # Step 2: convert to node class
    print("Generating demographics file and overlays")

    nodes = [Node(this_site["lat"], this_site["lon"],
          instructions["node_pop"],
          name= this_site["name"] if "name" in this_site.index else "",
          forced_id=this_site["cluster"],
          extra_attributes={"Country": this_site["birth_rate_country"]})
     for ix, this_site in sites.iterrows()]

    site_demog = DemographicsGenerator(nodes, res_in_arcsec="custom",
                               update_demographics=update_demog,
                               vectors= vector_props)

    demographics = site_demog.generate_demographics()
    demo_f = open(demog_out_path, "w+")
    json.dump(demographics, demo_f, indent=4)
    demo_f.close()

    overlay_path = os.path.join(demog_out_dir, "demographics_net_overlay.json")
    net_usage_overlay(demog_out_path, overlay_path)


