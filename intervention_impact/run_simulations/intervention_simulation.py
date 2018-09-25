import os
import pandas as pd
import numpy as np
import pdb
import json

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_asset_collection

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report, add_event_counter_report
from simtools.Utilities.Experiments import retrieve_experiment

from generate_input_files import generate_input_files, net_usage_overlay
from sweep_functions import *

# variables
run_type = "burnin"  # set to "burnin" or "intervention"
burnin_id = "8c066962-e6b5-e811-a2c0-c4346bcb7275"
asset_exp_id = "8c066962-e6b5-e811-a2c0-c4346bcb7275"
# intervention_coverages = [0, 20, 40, 60, 80]
# hs_daily_probs = [0.15, 0.3, 0.7]

int_scenarios = [{"act": 0,
                  "irs": 0,
                  "itn": 60},
                 {"act": 0,
                  "irs": 60,
                  "itn": 0},
                 {"act": 0,
                  "irs": 60,
                  "itn": 60},
                 {"act": 20,
                  "irs": 0,
                  "itn": 0},
                 {"act": 40,
                  "irs": 0,
                  "itn": 0},
                 {"act": 20,
                  "irs": 40,
                  "itn": 0},
                 {"act": 20,
                  "irs": 40,
                  "itn": 40},
                 ]

net_hating_props = [0.1] # based on expert opinion from Caitlin
new_inputs = True

# Serialization
print("setting up")
if run_type == "burnin":
    years = 15
    sweep_name = "MAP_II_Burnin_3"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    sweep_name = "MAP_II_Gates_Rerun"
    serialize = False
    pull_from_serialization = True
else:
    raise ValueError("Unknown run type " + run_type)

# setup
location = "HPC"
SetupParser.default_block = location


cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                    Simulation_Duration=int(365*years),
                                    Config_Name=sweep_name,
                                    Birth_Rate_Dependence="FIXED_BIRTH_RATE",
                                    Age_Initialization_Distribution_Type= "DISTRIBUTION_COMPLEX",
                                    Num_Cores=1,

                                    # interventions
                                    Valid_Intervention_States=[],  # apparently a necessary parameter
                                    # todo: do I need listed events?
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One", "Bednet_Using"],
                                    Enable_Default_Reporting=0,
                                    Enable_Demographics_Risk=1,

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
                  "Enable_Property_Output": 0})

if serialize:
    cb.update_params({"Serialization_Time_Steps": [365*years]})

# reporting
add_summary_report(cb)
# add_event_counter_report(cb, ["Bednet_Using"])

def set_asset_id(cb, asset_collection, use_assets=True):

    if use_assets:
        cb.set_input_collection(asset_collection)
        return {"Input_collection": str(asset_collection.id)}
    else:
        return {"Input_collection": "None"}

if __name__=="__main__":

    SetupParser.init()

    # collect site-specific data to pass to builder functions
    COMPS_login("https://comps.idmod.org")
    sites = pd.read_csv("site_details.csv")

    print("finding collection ids and vector details")
    # collection ids:
    cb.set_exe_collection("66483753-b884-e811-a2c0-c4346bcb7275")
    cb.set_dll_collection("65483753-b884-e811-a2c0-c4346bcb7275")

    site_input_dir = os.path.join("sites", "all")

    with open("species_details.json") as f:
        species_details = json.loads(f.read())

    if asset_exp_id:
        print("retrieving asset experiment")
        asset_expt = retrieve_experiment(asset_exp_id)
        asset_df = pd.DataFrame([x.tags for x in asset_expt.simulations])
        set_asset_id(cb, asset_expt)

    if new_inputs:
        print("generating input files")
        generate_input_files(site_input_dir, pop=2000, overwrite=True)

    # Find vector proportions for each vector in our site todo:use this for demog setup
    vectors = pd.read_csv(os.path.join(site_input_dir, "vector_proportions.csv"))

    simulation_setup(cb, species_details)

    if pull_from_serialization:
        print("building from pickup")

        # serialization
        print("retrieving burnin")
        expt = retrieve_experiment(burnin_id)

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        # temp for testing
        # df = df.query("Site_Name=='aba' &  Run_Number==7")

        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames": [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name],
                "Run_Number": df["Run_Number"][x],
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]}),

            ModFn(add_annual_itns, year_count=years,
                                   n_rounds=1,
                                   coverage= int_combo["itn"] / 100,
                                   discard_halflife=180,
                                   start_day=5,
                                   IP=[{"NetUsage":"LovesNets"}]
                  ),
            ModFn(assign_net_ip, hates_net_prop),
            ModFn(add_irs_group, coverage= int_combo["irs"]/100,
                                 decay=180,
                                 start_days=[365*start for start in range(years)]),
            ModFn(add_healthseeking_by_coverage, coverage=int_combo["act"]/100, rate=0.15),

        ]
            for x in df.index
            for hates_net_prop in net_hating_props
            for int_combo in int_scenarios
            # for daily_prob in hs_daily_probs

        ])

        run_sim_args = {"config_builder": cb,
                        "exp_name": sweep_name,
                        "exp_builder": builder}

        em = ExperimentManagerFactory.from_cb(cb)
        em.run_simulations(**run_sim_args)


    else:
        print("building burnin")
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Run_Number": run_num,
                "x_Temporary_Larval_Habitat":10 ** hab_exp}),
        ]
            for run_num in range(10)
            for hab_exp in np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))
        ])

    run_sim_args = {"config_builder": cb,
                    "exp_name": sweep_name,
                    "exp_builder": builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

