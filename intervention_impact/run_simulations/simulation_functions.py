###############################################################################################################
## simulation_functions.py
## Amelia Bertozzi-Villa
## October 2019
##
## Library of functions to set up intervention impact simulations.
##############################################################################################################

## Importing  ---------------------------------------------------------------------------------------

import pandas as pd
import json
import math

import os
import pdb

from simtools.Utilities.Experiments import retrieve_experiment

from dtk.vector.species import set_params_by_species, set_species_param
from dtk.interventions.habitat_scale import scale_larval_habitats
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.property_change import change_individual_property
from dtk.interventions.novel_vector_control import add_ATSB, add_ors_node, add_larvicides
from dtk.interventions.ivermectin import add_ivermectin

from malaria.reports.MalariaReport import add_summary_report
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from malaria.interventions.health_seeking import add_health_seeking

## Simulation setup ---------------------------------------------------------------------------------------

def assign_net_ip(cb, hates_net_prop):
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop),
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop,
          trigger_condition_list=["Births"])
    return {"Hates_Nets": hates_net_prop}


def set_up_simulation(cb, instructions, max_larval_capacity=4e8):

    # configure root directories for asset management
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cb.set_experiment_executable(os.path.join(dir_path, "bin", "Eradication.exe"))
    cb.set_dll_root(os.path.join(dir_path, "bin"))
    # cb.set_python_path()
    cb.set_input_files_root(instructions["root_dir"])

    # general settings
    print("customizing config")
    cb.update_params({"Birth_Rate_Depencence": "FIXED_BIRTH_RATE",
                      "Age_Initialization_Distribution_Type": "DISTRIBUTION_COMPLEX",

                      # interventions
                      "Listed_Events" : ["Bednet_Discarded", "Bednet_Got_New_One",
                                       "Bednet_Using", "Received_Vaccine"],
                      "Enable_Default_Reporting" : 0,
                      "Enable_Demographics_Risk" : 1,
                      "Enable_Vector_Species_Report": 0,

                      # ento from prashanth
                      "Antigen_Switch_Rate" : pow(10, -9.116590124),
                      "Base_Gametocyte_Production_Rate" : 0.06150582,
                      "Base_Gametocyte_Mosquito_Survival_Rate" : 0.002011099,
                      "Falciparum_MSP_Variants" : 32,
                      "Falciparum_Nonspecific_Types" :76,
                      "Falciparum_PfEMP1_Variants" : 1070,
                      "Gametocyte_Stage_Survival_Rate" : 0.588569307,
                      "MSP1_Merozoite_Kill_Fraction" : 0.511735322,
                      "Max_Individual_Infections" : 3,
                      "Nonspecific_Antigenicity_Factor" : 0.415111634,

                      # reporting
                      "Disable_IP_Whitelist": 1,
                      "Enable_Property_Output": 0

                      }
                     )

    # climate and demographic files
    demog_path = os.path.join(instructions["subdir_name"],  "demog")
    vector_path = os.path.join(instructions["subdir_name"], "vector")
    climate_path = os.path.join(instructions["subdir_name"],
                                "climate", instructions["this_run_type"])
    climate_path = climate_path if os.path.exists(os.path.join(instructions["root_dir"], climate_path)) else \
        os.path.join(climate_path, "..")
    climate_fnames = instructions["climate_fnames"]

    cb.update_params({
                    "Demographics_Filenames": [os.path.join(demog_path, "demographics.json"),
                                               os.path.join(demog_path, "demographics_net_overlay.json")],
                    "Air_Temperature_Filename": os.path.join(climate_path,
                                                             climate_fnames["air_temp"]),
                    "Land_Temperature_Filename": os.path.join(climate_path,
                                                              climate_fnames["land_temp"]),
                    "Rainfall_Filename": os.path.join(climate_path,
                                                      climate_fnames["rainfall"]),
                    "Relative_Humidity_Filename": os.path.join(climate_path,
                                                               climate_fnames["relative_humidity"])
                    }
    )

    # change net-hating proportion
    print("adding net ip")
    assign_net_ip(cb, 0.1)  # 10% based on opinion from Caitlin Bever

    # Find vector counts for each vector based on relative abundances
    print("adding vectors and scaling larval habitats")

    # custom list of species parameters
    species_details = { "arabiensis": {"habitat_split": {"CONSTANT": 0.1,
                                                    "TEMPORARY_RAINFALL": 0.9
                                                      },
                                    "Anthropophily": 0.65,
                                    "Indoor_Feeding_Fraction": 0.5
                                    },
                        "darlingi": {"habitat_split": {
                                                  "WATER_VEGETATION": 0.9,
                                                  "CONSTANT": 0.09,
                                                  "BRACKISH_SWAMP": 0.01
                                                 },
                                  "Anthropophily": 0.5
                                  },
                        "funestus": {"habitat_split": {
                                                    "WATER_VEGETATION": 1
                                                    },
                                  "Anthropophily": 0.65,
                                  "Indoor_Feeding_Fraction": 0.85

                                  },
                        "gambiae": {"habitat_split": {
                                                    "CONSTANT": 0.1,
                                                    "TEMPORARY_RAINFALL": 0.9
                                                   },
                                 "Anthropophily": 0.85,
                                 "Indoor_Feeding_Fraction": 0.85
                                 },
                        "maculatus": {"habitat_split": {
                                                    "TEMPORARY_RAINFALL": 0.8,
                                                    "WATER_VEGETATION": 0.1,
                                                    "CONSTANT": 0.1
                                                    },
                                   "Anthropophily": 0.5
                                   },
                        "minimus": {"habitat_split": {
                                                   "WATER_VEGETATION":  0.9,
                                                   "CONSTANT": 0.1
                                                  }
                        }
                    }

    site_vector_props = pd.read_csv(os.path.join(instructions["root_dir"], vector_path, "vector_proportions.csv"))
    larval_habs_per_site = {"NodeID": site_vector_props["id"]}
    species_needed = [value for value in species_details.keys() if value in site_vector_props.columns]
    set_params_by_species(cb.params, [name for name in species_needed])

    for species_name, species_modifications in species_details.items():
        if species_name not in species_needed:
            continue
        set_species_param(cb, species_name, "Adult_Life_Expectancy", 20)
        set_species_param(cb, species_name, "Vector_Sugar_Feeding_Frequency", "VECTOR_SUGAR_FEEDING_EVERY_DAY")

        for param, val in species_modifications.items():
            if param == "habitat_split":
                new_vals = {hab: hab_prop * max_larval_capacity for hab, hab_prop in val.items()}
                set_species_param(cb, species_name, "Larval_Habitat_Types", new_vals)
                larval_habs_per_site.update({".".join([species_name, hab]): site_vector_props[species_name]
                                             for hab in val.keys()})
            else:
                set_species_param(cb, species_name, param, val)

    # this function adds a per-node campaign event that scales larval habitat as specified.
    scale_larval_habitats(cb, pd.DataFrame(larval_habs_per_site))

    # reporting
    print("adding node-specific reports")
    # todo: can we use node demographics instead?
    for idx, row in site_vector_props.iterrows():
        add_summary_report(cb,
                           age_bins=list(range(10, 130, 10)),
                           nodes={
                               "class": "NodeSetNodeList",
                               "Node_List": [int(row["id"])]
                           },
                           description=str(row["id"]))

    # if applicable, set assets
    if instructions["asset_exp_id"]!="":
        print("retrieving asset experiment")
        asset_expt = retrieve_experiment(instructions["asset_exp_id"])
        template_asset = asset_expt.simulations[0].tags
        cb.set_exe_collection(template_asset["exe_collection_id"])
        cb.set_dll_collection(template_asset["dll_collection_id"])
        cb.set_input_collection(template_asset["input_collection_id"])

    print("simulation setup complete")




