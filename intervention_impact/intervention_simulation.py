import os
import pandas as pd
import numpy as np
import pdb

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.interventions.outbreakindividual import recurring_outbreak
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report
from sweep_functions import add_annual_itns, add_irs_group, add_healthseeking_by_coverage, site_simulation_setup

# variables
run_type = "burnin"  # set to "burnin" or "intervention"
# burnin_id = "d9101f31-1785-e811-a2c0-c4346bcb7275" # badly spaced burnin
burnin_id = "4683d052-b58a-e811-a2c0-c4346bcb7275"
intervention_coverages = [0, 20, 40, 60, 80]


# Serialization
if run_type == "burnin":
    years = 10
    exp_name = "MAP_II_Burnin_0"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    exp_name = "ACT_Impact_All_Sites"
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
                                    Age_Initialization_Distribution_Type= "DISTRIBUTION_COMPLEX",
                                    Num_Cores=1,

                                    # interventions
                                    Valid_Intervention_States=[],  # apparently a necessary parameter
                                    # todo: do I need listed events?
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One", "Bednet_Using"],
                                    Enable_Default_Reporting=0,

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

# collection ids:
cb.set_exe_collection("66483753-b884-e811-a2c0-c4346bcb7275")
cb.set_dll_collection("65483753-b884-e811-a2c0-c4346bcb7275")
def set_site_id(cb, site):
    input_collection_ids = {"aba": "75b33195-1685-e811-a2c0-c4346bcb7275",
                            "martae": "5db33195-1685-e811-a2c0-c4346bcb7275",
                            "djibo": "59b33195-1685-e811-a2c0-c4346bcb7275",
                            "kananga": "18aee64b-1685-e811-a2c0-c4346bcb7275",
                            "moine": "27b33195-1685-e811-a2c0-c4346bcb7275",
                            "gode": "b8eb358f-1685-e811-a2c0-c4346bcb7275",
                            "karen": "4feb358f-1685-e811-a2c0-c4346bcb7275",
                            "bajonapo": "49eb358f-1685-e811-a2c0-c4346bcb7275"}

    cb.set_input_collection(input_collection_ids[site])

# reporting
add_summary_report(cb)
cb.update_params({"Enable_Vector_Species_Report": 0})



if __name__=="__main__":
    SetupParser.init()

    if pull_from_serialization:

        # serialization
        COMPS_login("https://comps.idmod.org")
        expt = ExperimentDataStore.get_most_recent_experiment(burnin_id)

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        # temp for testing
       #  df = df.query("Site_Name=='moine'")

        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path", os.path.join(df["outpath"][x], "output")),
            ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Filenames",
                  [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name]),
            ModFn(DTKConfigBuilder.set_param, "Run_Number", df["Run_Number"][x]),
            ModFn(DTKConfigBuilder.set_param, "x_Temporary_Larval_Habitat", df["x_Temporary_Larval_Habitat"][x]),
            ModFn(set_site_id, site=df["Site_Name"][x]),
            ModFn(site_simulation_setup, site_name=df["Site_Name"][x]),

            # ModFn(add_annual_itns, year_count=years, n_rounds=1, coverage=itn_cov / 100, discard_halflife=180),
            # ModFn(recurring_outbreak, outbreak_fraction=outbreak_fraction, repetitions=12 * years, tsteps_btwn=30),
            # ModFn(add_irs_group, coverage=irs_cov/100, decay=180, start_days=[365*start for start in range(years)]),
            ModFn(add_healthseeking_by_coverage, coverage=act_cov/100),

        ]
            for x in df.index
            # for itn_cov in intervention_coverages
            # for n_dists in [1,2,3]
            # for outbreak_fraction in [0.001, 0.005, 0.01]
            # for irs_cov in intervention_coverages
            for act_cov in intervention_coverages

        ])
    else:
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.set_param, "x_Temporary_Larval_Habitat", 10 ** x),
            ModFn(DTKConfigBuilder.set_param, "Run_Number", y),
            ModFn(site_simulation_setup, site_name=site_name),
            ModFn(set_site_id, site=site_name)
        ]
            for x in np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))
            for y in range(10)
            for site_name in sites["name"]
        ])

    run_sim_args = {"config_builder": cb,
                    "exp_name": exp_name,
                    "exp_builder": builder}


    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

