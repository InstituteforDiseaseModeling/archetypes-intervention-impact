import os
import pandas as pd
import numpy as np
import pdb
import json
import math

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_asset_collection

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report, add_event_counter_report
from dtk.utils.reports.VectorReport import add_vector_stats_report
from simtools.Utilities.Experiments import retrieve_experiment

# from generate_input_files import generate_input_files, net_usage_overlay
from malaria.interventions.malaria_vaccine import add_vaccine
from sweep_functions import *

# variables
run_type = "intervention"  # set to "burnin" or "intervention"
burnin_id = "96e9c858-a8ce-e811-a2bd-c4346bcb1555"
asset_exp_id = "66d8416c-9fce-e811-a2bd-c4346bcb1555"

intervention_coverages = [60]
vaccine_durations = [182, 365]
interventions = ["itn_irs", "itn", 'irs', 'itn_w_irs', 'al_cm']
# hs_daily_probs = [0.15, 0.3, 0.7]

net_hating_props = [0.1] # based on expert opinion from Caitlin
new_inputs = False

# Serialization
print("setting up")
if run_type == "burnin":
    years = 15
    sweep_name = "MAP_II_New_Sites_Burnin"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    sweep_name = "ITN_IRS_with_cases"
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
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One", "Bednet_Using", "Received_Vaccine"],
                                    Enable_Default_Reporting=1,
                                    Enable_Demographics_Risk=1,
                                    Enable_Vector_Species_Report=0,

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
                  "Enable_Property_Output": 0,
                  "Enable_Spatial_Output": 1,
                  "Spatial_Output_Channels": ["Population", "Blood_Smear_Parasite_Prevalence", 'New_Infections',
                                              'New_Clinical_Cases']
                  })

if serialize:
    cb.update_params({"Serialization_Time_Steps": [365*years]})


if __name__=="__main__":

    SetupParser.init()

    # collect site-specific data to pass to builder functions
    COMPS_login("https://comps.idmod.org")
    sites = pd.read_csv("site_details.csv")

    print("finding collection ids and vector details")
    site_input_dir = os.path.join("sites", "all")

    with open("species_details.json") as f:
        species_details = json.loads(f.read())

    if asset_exp_id:
        print("retrieving asset experiment")
        asset_expt = retrieve_experiment(asset_exp_id)
        template_asset = asset_expt.simulations[0].tags
        cb.set_exe_collection(template_asset["exe_collection_id"])
        cb.set_dll_collection(template_asset["dll_collection_id"])
        cb.set_input_collection(template_asset["input_collection_id"])

    # if new_inputs:
    #     print("generating input files")
    #     generate_input_files(site_input_dir, pop=2000, overwrite=True)

    # Find vector proportions for each vector in our site
    site_vectors = pd.read_csv(os.path.join(site_input_dir, "vector_proportions.csv"))
    simulation_setup(cb, species_details, site_vectors)

    # reporting
    for idx, row in site_vectors.iterrows():
        add_summary_report(cb,
                           age_bins = list(range(10, 130, 10)),
                           nodes={
                               "class": "NodeSetNodeList",
                               "Node_List": [int(row["node_id"])]
                           },
                           description = row["name"])
    # add_event_counter_report(cb, ["Bednet_Using", "Received_Vaccine"])
    # add_vector_stats_report(cb)

    if pull_from_serialization:
        print("building from pickup")

        # serialization
        print("retrieving burnin")
        expt = retrieve_experiment(burnin_id)

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        # temp for testing
        # df = df.query("Run_Number==7 & x_Temporary_Larval_Habitat>90 & x_Temporary_Larval_Habitat<150")

        old_builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames": [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name],
                "Run_Number": df["Run_Number"][x],
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]}),

            ModFn(add_annual_itns, year_count=years,
                                   n_rounds=1,
                                   coverage= itn_cov / 100,
                                   discard_halflife=180,
                                   start_day=5,
                                   IP=[{"NetUsage":"LovesNets"}]
                  ),
            ModFn(assign_net_ip, hates_net_prop),
            ModFn(add_irs_group, coverage= irs_cov /100,
                                 decay=180,
                                 start_days=[365*start for start in range(years)]),
            ModFn(add_healthseeking_by_coverage, coverage= act_cov/100, rate=0.15),

        ]
            for x in df.index
            for hates_net_prop in net_hating_props
            for itn_cov in intervention_coverages
            for irs_cov in intervention_coverages
            for act_cov in intervention_coverages

        ])

        from_burnin_list = [
            ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames": [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name],
                "Run_Number": df["Run_Number"][x],
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]})

        for x in df.index]

        intervention_dict = {
            "itn": list(
            [burnin_fn,
             ModFn(add_annual_itns, year_count=years,
                   n_rounds=1,
                   coverage=itn_cov / 100,
                   discard_halflife=180,
                   start_day=5,
                   IP=[{"NetUsage": "LovesNets"}]
                   ),
             ModFn(assign_net_ip, hates_net_prop),
             ModFn(add_healthseeking_by_coverage, coverage=40 / 100, rate=0.15)
             ]

            for burnin_fn in from_burnin_list
            for itn_cov in intervention_coverages
            for hates_net_prop in net_hating_props
            ),

            "itn_irs": list(
                [burnin_fn,
                 ModFn(add_annual_itns, year_count=years,
                       n_rounds=1,
                       coverage=itn_cov / 100,
                       discard_halflife=180,
                       start_day=5,
                       IP=[{"NetUsage": "LovesNets"}]
                       ),
                 ModFn(assign_net_ip, hates_net_prop),
                 ModFn(add_irs_group, coverage=itn_cov / 100,
                       decay=180,
                       start_days=[365 * start for start in range(years)]),
                 ModFn(add_healthseeking_by_coverage, coverage=40 / 100, rate=0.15)
                 ]

                for burnin_fn in from_burnin_list
                for itn_cov in intervention_coverages
                for hates_net_prop in net_hating_props
            ),

            "irs": list(
            [burnin_fn,
             ModFn(add_irs_group, coverage=irs_cov / 100,
                   decay=180,
                   start_days=[365 * start for start in range(years)]),
             ModFn(add_healthseeking_by_coverage, coverage=40 / 100, rate=0.15)
             ]
            for burnin_fn in from_burnin_list
            for irs_cov in intervention_coverages
            ),

            "itn_w_irs": list(
                [burnin_fn,
                 ModFn(add_annual_itns_w_irs, year_count=years,
                       n_rounds=1,
                       coverage=itn_cov / 100,
                       discard_halflife=180,
                       start_day=5,
                       IP=[{"NetUsage": "LovesNets"}]
                       ),
                 ModFn(add_healthseeking_by_coverage, coverage=40 / 100, rate=0.15)
                 ]
                for burnin_fn in from_burnin_list
                for itn_cov in intervention_coverages
            ),

            "al_cm": list(
            [burnin_fn,
             ModFn(add_healthseeking_by_coverage, coverage=al_cm_cov / 100, rate=0.15)
             ]
            for burnin_fn in from_burnin_list
            for al_cm_cov in [40]#intervention_coverages
            ),

            "dp_cm": list(
                [burnin_fn,
                 ModFn(add_healthseeking_by_coverage, coverage=dp_cm_cov / 100, rate=0.15, drugname="DP")
                 ]
                for burnin_fn in from_burnin_list
                for dp_cm_cov in intervention_coverages
            ),

            "dp_mda": list(
                [burnin_fn,
                 ModFn(add_mda, coverage=dp_mda_cov/100)]

                for burnin_fn in from_burnin_list
                for dp_mda_cov in intervention_coverages
            ),

            "mAb":list(
                [burnin_fn,
                 ModFn(add_vaccine, vaccine_type="PEV",
                                    coverage=mab_cov/100,
                                    vaccine_params={"Waning_Config": {
                                                        "class": "WaningEffectBox",
                                                        "Box_Duration": 90 # 3 month protection
                                                    }
                                                    },
                                    target_group={"agemin": 15, "agemax": 49},
                                    repetitions=3,
                                    interval=365
                       )
                ]

                for burnin_fn in from_burnin_list
                for mab_cov in intervention_coverages
            ),

            "pev": list(
                [burnin_fn,
                 ModFn(add_vaccine, vaccine_type="PEV",
                                    coverage=pev_cov/100,
                                    vaccine_params={"Reduced_Acquire": 0.75,
                                                    "Waning_Config": {
                                                        "class": "WaningEffectExponential",
                                                        "Decay_Time_Constant": vaccine_hl/ math.log(2)
                                                    }
                                                    },
                                    trigger_condition_list=["Births"],
                                    triggered_delay=182
                       )
                ]

                for burnin_fn in from_burnin_list
                for pev_cov in intervention_coverages
                for vaccine_hl in vaccine_durations
            ),

            "tbv": list(
                [burnin_fn,
                 ModFn(add_vaccine, vaccine_type="TBV",
                       coverage=tbv_cov/100,
                       vaccine_params={"Reduced_Transmit": 0.75,
                                       "Waning_Config": {
                                           "class": "WaningEffectExponential",
                                           "Decay_Time_Constant": vaccine_hl / math.log(2)
                                       }
                                       },
                       target_group={"agemin": 15, "agemax": 49},
                       repetitions = 3,
                       interval = 365
                       )
                 ]

                for burnin_fn in from_burnin_list
                for tbv_cov in intervention_coverages
                for vaccine_hl in vaccine_durations
            )

        }

        combos = []
        for int in interventions:
            combos.extend(intervention_dict[int])

        builder = ModBuilder.from_list(combos)


    else:
        print("building burnin")
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Run_Number": run_num,
                "x_Temporary_Larval_Habitat":10 ** hab_exp}),
        ]
            for run_num in range(10)
            for hab_exp in np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))
            # for hab_exp in [0, 1, 2]
        ])

    run_sim_args = {"config_builder": cb,
                    "exp_name": sweep_name,
                    "exp_builder": builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

