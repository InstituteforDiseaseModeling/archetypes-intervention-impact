###############################################################################################################
## 02_generate_climate.py
## Amelia Bertozzi-Villa
## October 2019
##
## EMOD needs climate files and demographics files. Based on the instructions specified in the input_params.json,
## this script creates ERA5 climate parameters for the specified times and locations.

## Requirements: input directory containing an "input_params.json" and a "site_details.csv" file.
##                see "example_input_dir" in this repo.
##
## Outputs: A directory named "climate", possibly with subdirectories "burnin" and "intervention"
## (if different climate files are needed for each). Relevant input files must be collected from COMPS
## and stored in these directories.
##############################################################################################################

## Importing and setup ---------------------------------------------------------------------------------------
import pandas as pd
import os
import json

from simtools.AssetManager.FileList import FileList
from simtools.Managers.WorkItemManager import WorkItemManager
from simtools.SetupParser import SetupParser

desired_width = 320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)

## Function ---------------------------------------------------------------------------------------------

def check_files(fname_dict, base_dir):
    files_exist = [os.path.exists(os.path.join(base_dir, this_fname))
                    for this_fname in fname_dict.values()]
    return(all(files_exist))


def confirm_climate_files(instructions, climate_dir):
    if "era5_climate_params" in instructions.keys():
        all_files_exist = [check_files(instructions["climate_fnames"],
                                       os.path.join(climate_dir, subdir))
                           for subdir in instructions["era5_climate_params"].keys()]
    else:
        all_files_exist = [check_files(instructions["climate_fnames"], climate_dir)]
    return(all(all_files_exist))

## Main code setup ---------------------------------------------------------------------------------------------

main_dir = os.path.join(os.path.expanduser("~"),
                        "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                        "20191009_megatrends_era5_new_archetypes",
                        "input")

# uncomment this to run the example
# main_dir = "example_input_dir"

with open(os.path.join(main_dir, "input_params.json")) as f:
    instructions = json.loads(f.read())

main_input_dir = instructions["root_dir"]

# A .csv or a demographics file containing input coordinates
site_fname = instructions["site_fname"]
sites = pd.read_csv(os.path.join(main_dir, site_fname))

## ERA5 Climate generation ---------------------------------------------------------------------------------------
# (or, if other climate files are already generated, do nothing).

climate_out_dir = os.path.join(main_dir, "climate")
climate_files_exist = confirm_climate_files(instructions, climate_out_dir)

if climate_files_exist and instructions["overwrite_input_files"] == "False":
    print("Climate files already generated")
else:
    SetupParser.default_block = 'HPC'
    SetupParser.init()

    for run_type, years in instructions["era5_climate_params"].items():
        print("Generating climate for {run_type},  {start_year} to {end_year}".format(run_type=run_type,
                                                                                      start_year=years["start_year"],
                                                                                      end_year=years["end_year"]))
        this_out_dir = os.path.join(climate_out_dir, run_type)
        if not os.path.isdir(this_out_dir):
            os.makedirs(this_out_dir)

        wi_name = "ERA5 weather: {run_type}".format(run_type=run_type)
        optional_args = "--ds ERA5 --id-ref 'Custom user' --node-col cluster"

        # To run a specific version add a tag (for example, "weather-files:1.1").
        # See available versions here: https://github.com/InstituteforDiseaseModeling/dst-era5-weather-data-tools/releases.
        docker_image = "weather-files"
        command_pattern = "python /app/generate_weather_asset_collection.py {} {} {} {}"
        command = command_pattern.format(site_fname, years["start_year"], years["end_year"], optional_args)
        user_files = FileList(root=main_dir, files_in_root=[site_fname])

        wi = WorkItemManager(item_name=wi_name, docker_image=docker_image, command=command, user_files=user_files,
            tags={"Demo": "dtk-tools Docker WorkItem", "WorkItem type": "Docker", "Command": command })

        wi.create()
        wi.run()

    print("Climate file generation submitted. Go to COMPS to monitor and retrieve them.")