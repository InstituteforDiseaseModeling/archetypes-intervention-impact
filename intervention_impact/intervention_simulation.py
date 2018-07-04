import os
import pandas as pd
import numpy as np

from simtools.SetupParser import SetupParser

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report
from sweep_functions import add_annual_itns, add_irs_group, add_healthseeking_by_coverage, site_simulation_setup

# variables
run_type = "burnin"  # set to "burnin" or "intervention"
burnin_id = "b5f4c942-b67f-e811-a2c0-c4346bcb7275"

# Serialization
if run_type == "burnin":
    years = 50
    exp_name = "Intervention_Impact_Burnins"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    exp_name = "Intervention_Impact_Interventions"
    serialize = False
    pull_from_serialization = True
else:
    raise ValueError("Unknown run type " + run_type)

# setup
location = "HPC"
SetupParser.default_block = location
sites = pd.read_csv("site_details.csv")

cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                    Simulation_Duration=int(365*years),
                                    Config_Name=exp_name,
                                    Birth_Rate_Dependence="FIXED_BIRTH_RATE",
                                    Num_Cores=1,

                                    # interventions
                                    Valid_Intervention_States=[],  # apparently a necessary parameter
                                    # todo: do I need listed events?
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One", "Bednet_Using"],

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

if serialize:
    cb.update_params({"Serialization_Time_Steps": [365*years]})

# reporting
add_summary_report(cb)

if pull_from_serialization:

    # serialization
    COMPS_login("https://comps.idmod.org")
    expt = ExperimentDataStore.get_most_recent_experiment(burnin_id)

    df = pd.DataFrame([x.tags for x in expt.simulations])
    df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path", os.path.join(df["outpath"][x], "output")),
        ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Filenames",
              [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name]),
        ModFn(DTKConfigBuilder.set_param, "Run_Number", df["Run_Number"][x]),
        ModFn(DTKConfigBuilder.set_param, "x_Temporary_Larval_Habitat", df["x_Temporary_Larval_Habitat"][x]),
        ModFn(add_annual_itns, year_count=years, n_rounds=1, coverage=itn_cov/100, discard_halflife=180),
        ModFn(add_healthseeking_by_coverage, coverage=irs_cov/100),
        ModFn(add_irs_group, coverage=act_cov/100, decay=180, start_days=[365*start for start in range(years)]),
        ModFn(site_simulation_setup, site_name=df["Site_Name"][x])
                                     ]
        for x in df.index
        for itn_cov in [0, 80]
        for irs_cov in [0, 80]
        for act_cov in [0, 80]

    ])
else:
    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.set_param, "x_Temporary_Larval_Habitat", 10**x),
        ModFn(DTKConfigBuilder.set_param, "Run_Number", y),
        ModFn(site_simulation_setup, site_name=site_name)
        ]
        for x in np.concatenate((np.arange(0, 2.25, 0.05), np.arange(2.25, 4.25, 0.25)))
        for y in range(5)
        for site_name in sites["name"]
    ])

run_sim_args = {"config_builder": cb,
                "exp_name": exp_name,
                "exp_builder": builder}
