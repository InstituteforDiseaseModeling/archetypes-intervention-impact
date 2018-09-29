import os
import numpy as np

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_asset_collection

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.vector.species import set_params_by_species, set_species_param
from simtools.ModBuilder import ModBuilder, ModFn

from malaria.reports.MalariaReport import add_summary_report, add_event_counter_report


# setup-- these variables change often. Unless you change the "builder" parameters below,
# "years" should be at least 2!
years = 2
exp_name = "example_run"
location = "HPC"
SetupParser.default_block = location

# set up a basic config with the latest entomology params
cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                    Simulation_Duration=int(365*years),
                                    Config_Name=exp_name,
                                    Birth_Rate_Dependence="FIXED_BIRTH_RATE",
                                    Age_Initialization_Distribution_Type= "DISTRIBUTION_COMPLEX",
                                    Num_Cores=1,

                                    # interventions
                                    Valid_Intervention_States=[],  # apparently a necessary parameter
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One", "Bednet_Using"],
                                    Enable_Default_Reporting=1,
                                    Enable_Vector_Species_Report=1,

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

# update vector parameters
vectors = ["arabiensis", "funestus", "gambiae"]
set_params_by_species(cb.params, vectors)

for species_name in vectors:
    set_species_param(cb, species_name, "Adult_Life_Expectancy", 20)
    set_species_param(cb, species_name, "Vector_Sugar_Feeding_Frequency", "VECTOR_SUGAR_FEEDING_NONE")

# update input file paths (adjust as needed for your setup)
site_name = "aba"
site_dir = "sites/aba"
cb.update_params({
                "Demographics_Filenames": ["sites/{site}/demographics_{site}.json".format(site=site_name)],
                "Air_Temperature_Filename": os.path.join(site_dir,
                                                         "air_temperature_daily.bin"),
                "Land_Temperature_Filename": os.path.join(site_dir,
                                                          "air_temperature_daily.bin"),
                "Rainfall_Filename": os.path.join(site_dir,
                                                  "rainfall_daily.bin"),
                "Relative_Humidity_Filename": os.path.join(site_dir,
                                                           "relative_humidity_daily.bin")
                }
)

# add custom reporting
add_summary_report(cb, age_bins=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120])
add_event_counter_report(cb, ["Bednet_Using"])

# collection ids give a lot of speedup if you're using the same inputs over and over again--
# they pull .exes, .dlls, or input files from other simulations so you don't have to re-upload every time.
# I've commented these out b/c I don't know that you'll be able to access my experiments, but try them out!
# cb.set_exe_collection("66483753-b884-e811-a2c0-c4346bcb7275")
# cb.set_dll_collection("65483753-b884-e811-a2c0-c4346bcb7275")
# cb.set_input_collection("51448321-0291-e811-a2c0-c4346bcb7275")


# creates a bednet intervention and returns a dictionary that will get passed to COMPS as "tags"
def add_annual_itns(cb, year_count=1, coverage=0.8, discard_halflife=270, start_day=0, IP=[]):

    for year in range(year_count):

        add_ITN_age_season(cb,
                           coverage_all=coverage,
                           discard={"halflife": discard_halflife},
                           start=(365 * year) + start_day,
                           ind_property_restrictions=IP)

    return {"ITN_Coverage": coverage, "ITN_Halflife": discard_halflife,
            "ITN_start": start_day, "ITN_Distributions": year_count}


### SLOW WAY: ---------------------------------------------------------------------------------------------------------
# # The slow way, set up your simulation as below and run using "dtk run sample_sim_run.py"
#
# # Here's where you define what will make one simulation different from another. In this case I run:
# # 10 random seeds
# # 20 initial larval habitat magnitudes
# # 4 intervention coverages
# # for a total of 800 simulations.
#
# # this simulation will let things run for two years, then deploy bednets annually at the specified coverage
# # for the remainder of the run.
# builder = ModBuilder.from_list([[
#     ModFn(DTKConfigBuilder.update_params, {
#         "Run_Number": run_num,
#         "x_Temporary_Larval_Habitat":10 ** hab_exp}),
#     ModFn(add_annual_itns, year_count=years-2,
#           coverage=itn_cov / 100,
#           discard_halflife=180,
#           start_day=760
#           ),
# ]
#     for run_num in range(10)
#     for hab_exp in np.arange(0, 2, 0.1)
#     for itn_cov in [0, 20, 60, 80]
# ])
#
# # finally, this pulls everything together for "dtk run" to pick up
# run_sim_args = {"config_builder": cb,
#                 "exp_name": exp_name,
#                 "exp_builder": builder}
#

### --------------------------------------------------------------------------------------------------------------------

### BETTER WAY: ---------------------------------------------------------------------------------------------------------
# because of how dtk-tools interacts with COMPS, it is DRAMATICALLY FASTER to set up your script like
# this instead and simply call it using "python sample_sim_run.py" :

if __name__=="__main__":

    SetupParser.init()

    # Here's where you define what will make one simulation different from another. In this case I run:
    # 10 random seeds
    # 20 initial larval habitat magnitudes
    # 4 intervention coverages
    # for a total of 800 simulations.

    # this simulation will let things run for two years, then deploy bednets annually at the specified coverage
    # for the remainder of the run.
    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.update_params, {
            "Run_Number": run_num,
            "x_Temporary_Larval_Habitat":10 ** hab_exp}),
        ModFn(add_annual_itns, year_count=years-2,
              coverage=itn_cov / 100,
              discard_halflife=180,
              start_day=760
              ),
    ]
        for run_num in range(1)
        for hab_exp in [0,1,2]
        for itn_cov in [20, 80]
    ])

    # finally, this pulls everything together for "dtk run" to pick up
    run_sim_args = {"config_builder": cb,
                    "exp_name": exp_name,
                    "exp_builder": builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

### --------------------------------------------------------------------------------------------------------------------
