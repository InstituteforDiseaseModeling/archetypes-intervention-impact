import os
import pandas as pd
import numpy as np
import pdb
import json

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_asset_collection

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.interventions.outbreakindividual import recurring_outbreak
from dtk.interventions.biting_risk import change_biting_risk
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report, add_event_counter_report
# todo: fix fiona?? bug
# from generate_input_files import generate_input_files, net_usage_overlay
from sweep_functions import *

# variables
run_type = "intervention"  # set to "burnin" or "intervention"
# burnin_id = "d9101f31-1785-e811-a2c0-c4346bcb7275" # badly spaced burnin
burnin_id = "5553b581-d18f-e811-a2c0-c4346bcb7275"
intervention_coverages = [0, 20, 40, 60, 80]
intervention_coverages = [80]
net_hating_props = [0, 0.2, 0.5, 0.8]
new_inputs = False

# Serialization
if run_type == "burnin":
    years = 10
    exp_name = "MAP_II_Burnin_Test_IP"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    exp_name = "Test_Biting_Risk"
    serialize = False
    pull_from_serialization = True
else:
    raise ValueError("Unknown run type " + run_type)

# setup
location = "HPC"
SetupParser.default_block = location


cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                    Simulation_Duration=int(365*years),
                                    Config_Name=exp_name,
                                    Birth_Rate_Dependence="FIXED_BIRTH_RATE",
                                    Age_Initialization_Distribution_Type= "DISTRIBUTION_COMPLEX",
                                    Num_Cores=1,

                                    # interventions
                                    Valid_Intervention_States=[],  # apparently a necessary parameter
                                    # todo: do I need listed events?
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One", "Bednet_Using"],
                                    Enable_Default_Reporting=1,

                                    # ento from prashanth
                                    Antigen_Switch_Rate=pow(10, -9.116590124),
                                    Base_Gametocyte_Production_Rate=0.06150582,
                                    Base_Gametocyte_Mosquito_Survival_Rate=0.002011099,
                                    Falciparum_MSP_Variants=32,
                                    Falciparum_Nonspecific_Types=76,
                                    Falciparum_PfEMP1_Variants=1070,
                                    Gametocyte_Stage_Survival_Rate=0.588569307,
                                    MSP1_Merozoite_Kill_Fraction=0.511735322,
                                    Max_Individual_Infections=3,
                                    Nonspecific_Antigenicity_Factor=0.415111634,

                                    )

cb.update_params({"Disable_IP_Whitelist": 1,
                  "Enable_Property_Output": 1})

# add hetero biting
change_biting_risk(cb, risk_config={'Risk_Distribution_Type': 'EXPONENTIAL_DURATION', 'Exponential_Mean': 1})

if serialize:
    cb.update_params({"Serialization_Time_Steps": [365*years]})

# reporting
add_summary_report(cb)
add_event_counter_report(cb, ["Bednet_Using"])

def set_site_id(cb, asset_collection):
    cb.set_input_collection(asset_collection)
    return {"Input_collection": str(asset_collection.id)}

if __name__=="__main__":

    # collect site-specific data to pass to builder functions
    COMPS_login("https://comps.idmod.org")
    sites = pd.read_csv("site_details.csv")

    # collection ids:
    cb.set_exe_collection("66483753-b884-e811-a2c0-c4346bcb7275")
    cb.set_dll_collection("17f8bb9c-6f8f-e811-a2c0-c4346bcb7275")

    site_info = {}

    with open("species_details.json") as f:
        species_details = json.loads(f.read())

    for site_name in sites["name"]:
        site_dir = os.path.join("sites", site_name)
        site_info[site_name] = {}

        # input files
        if new_inputs:
            print("generating input files for " + site_name)
            # generate_input_files(site_name, pop=2000, overwrite=True)

        # make sure net overlay exists
        overlay_fname = "demographics_{name}_hatenets_0.json".format(name=site_name)
        # if not os.path.isfile(os.path.join("sites", site_name, overlay_fname)):
        #     net_usage_overlay(site_name, hates_net_prop=0)

        # asset collections
        site_info[site_name]["asset_collection"] = get_asset_collection(
            sites.query('name==@site_name')['asset_id'].iloc[0])

        # Find vector proportions for each vector in our site
        vectors = pd.read_csv(os.path.join(site_dir, "vector_proportions.csv"))
        site_info[site_name]["vectors"] = {row.species: row.proportion for row in vectors.itertuples() if
                                           row.proportion > 0}

    # builders
    SetupParser.init()

    if pull_from_serialization:

        # serialization
        expt = ExperimentDataStore.get_most_recent_experiment(burnin_id)

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        # temp for testing
        df = df.query("Site_Name=='moine' & x_Temporary_Larval_Habitat>10 & x_Temporary_Larval_Habitat<11 & Run_Number==9")

        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames": [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name],
                "Run_Number": df["Run_Number"][x],
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]}),
            # ModFn(set_site_id, asset_collection=site_info[df["Site_Name"][x]]["asset_collection"]),
            ModFn(site_simulation_setup, site_name=df["Site_Name"][x],
                                         species_details=species_details,
                                         vectors=site_info[df["Site_Name"][x]]["vectors"]),

            ModFn(add_annual_itns, year_count=years,
                                   n_rounds=1,
                                   coverage=itn_cov / 100,
                                   discard_halflife=180,
                                   start_day=5,
                                   IP=[{"NetUsage":"LovesNets"}]
                  ),
            ModFn(assign_net_ip, hates_net_prop),
            # ModFn(recurring_outbreak, outbreak_fraction=outbreak_fraction,
            #                           repetitions=12 * years,
            #                           tsteps_btwn=30),
            # ModFn(add_irs_group, coverage=irs_cov/100,
            #                      decay=180,
            #                      start_days=[365*start for start in range(years)]),
            # ModFn(add_healthseeking_by_coverage, coverage=act_cov/100),

        ]
            for x in df.index
            for itn_cov in intervention_coverages
            for hates_net_prop in net_hating_props
            # for n_dists in [1,2,3]
            # for outbreak_fraction in [0.001, 0.005, 0.01]
            # for irs_cov in intervention_coverages
            # for act_cov in intervention_coverages

        ])
    else:
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Run_Number": run_num,
                "x_Temporary_Larval_Habitat":10 ** hab_exp}),
            ModFn(set_site_id, asset_collection=site_info[site_name]["asset_collection"]),
            ModFn(site_simulation_setup, site_name=site_name,
                                         species_details=species_details,
                                         vectors=site_info[site_name]["vectors"]),
        ]
            for run_num in range(10)
            for hab_exp in np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))
            for site_name in sites["name"]
        ])

    run_sim_args = {"config_builder": cb,
                    "exp_name": exp_name,
                    "exp_builder": builder}


    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

