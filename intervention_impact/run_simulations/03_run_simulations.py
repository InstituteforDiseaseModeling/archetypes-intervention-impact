###############################################################################################################
## 03_run_simulations.py
## Amelia Bertozzi-Villa
## October 2019
##
## Using the climate, demographic, and (possibly) intervention inputs prepared in steps 0-2, this script submits
## the appropriate burnin or intervention simulations, pulling further input details from the "input_params.json"
## file in the input directory.
##
## Requirements: input directory containing:
##                       - an "input_params.json"
##                       - a "site_details.csv"
##                       - "vector", "demog", and "climate" directories with the appropriate inputs.
##                       - if running interventions, an "interventions.csv"
##  For formatting, see "example_input_dir" in this repo.

## Outputs: Submitted experiments with the expected behavior.
##############################################################################################################

## Importing and setup ---------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import itertools
import pdb
import os
import json

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.Utilities.Experiments import retrieve_experiment
from simtools.Utilities.COMPSUtilities import COMPS_login

from simulation_functions import *
from intervention_functions import *

vars = sorted(os.environ.keys())
print([os.environ[var] for var in vars])

desired_width = 320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)

## Workaround for strange asset issues on mac
os.environ['NO_PROXY'] = 'comps.idmod.org'

## VARIABLES-- user should set these ---------------------------------------------------------------------------------

node_group = "emod_abcd"
priority = "Lowest"

version_name = "20191009_megatrends_era5_new_archetypes"
main_dir = os.path.join(os.path.expanduser("~"),
                            "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                            version_name)
experiment_root_name = "MAP_" + version_name

# uncomment this to run the example
# main_dir = "example_input_dir"
# experiment_root_name = "int_impact_example"

# run_type: set to "burnin" or "intervention".
# If "intervention", the "burnin_id" field of "input_params.json" must be populated.
run_type = "burnin"
test_run = False
node_group = "emod_32cores" if test_run else "emod_abcd"

## Main code setup ---------------------------------------------------------------------------------------------------

# setup
location = "HPC"
SetupParser.default_block = location

if __name__=="__main__":
    SetupParser.init()
    SetupParser.set("HPC", "priority", priority)
    SetupParser.set("HPC", "node_group", node_group)

    COMPS_login("https://comps.idmod.org")

    with open(os.path.join(main_dir, "input_params.json")) as f:
        instructions = json.loads(f.read())

    sites = pd.read_csv(os.path.join(main_dir, instructions["site_fname"]))
    instructions["this_run_type"] = run_type

    # Serialization and naming
    if run_type == "burnin":
        years = instructions["burnin_years"]
        experiment_name = experiment_root_name + "_Burnin"
    elif run_type == "intervention":
        years = instructions["intervention_years"]
        experiment_name = experiment_root_name + "_Intervention"
    else:
        raise ValueError("Unknown run type " + run_type)

    experiment_name = "{experiment_name}_TEST".format(experiment_name=experiment_name) if test_run else experiment_name

    # initialize cb
    cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                        Simulation_Duration=int(365 * years),
                                        Config_Name=experiment_name
                                        )
    # run main setup function
    set_up_simulation(cb, instructions)

    ## Set up burnin ---------------------------------------------------------------------------------------------
    if run_type=="burnin":

        print("Building burnin")
        if test_run:
            print("Running three test sims")
            run_count = 1
            hab_exps = [0, 1, 2]
        else:
            run_count = 10
            hab_exps = np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))

        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Run_Number": run_num,
                "x_Temporary_Larval_Habitat": 10 ** hab_exp,
                "Serialization_Time_Steps": [365 * years]
            }),
        ]
            for run_num in range(run_count)
            for hab_exp in hab_exps
        ])

    ## Set up intervention scenarios --------------------------------------------------------------------------------
    else:

        print("Building intervention scenarios")

        # serialization
        print("retrieving burnin")
        expt = retrieve_experiment(instructions["burnin_id"])

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        if test_run:
            print("Running three test sims")
            df = df.iloc[0:1]

        # find burnin length for filename (should be the same for all sims in df)
        burnin_length_in_days = df["Serialization_Time_Steps"][0].strip('[]')

        from_burnin_list = [
            [ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames":
                    ["state-{burnin_str}.dtk".format(burnin_str=str(burnin_length_in_days).zfill(5))],
                "Run_Number": df["Run_Number"][x],
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]})]
            for x in df.index]

        # generate interventions
        def get_combos_and_flatten(old_list):
            combos = list(itertools.product(*old_list))
            return [list(itertools.chain.from_iterable(sublist)) for sublist in combos]

        # load intervention dataset
        interventions = pd.read_csv(os.path.join(main_dir, "interventions.csv"))
        intervention_dict = generate_intervention_tuples(coverages=interventions["cov"].unique().tolist(),
                                                         start_days=interventions["start_day"].unique().tolist(),
                                                         years=years)

        # generate a list of intervention packages from dataset
        # (each sublist is its own intervention package)
        full_int_list = []
        for this_int_idx in interventions["int_id"].unique().tolist():
            this_int_package = interventions.query("int_id==@this_int_idx")
            this_int_list = [intervention_dict[row["start_day"]][row["cov"]][row["int"]]
            for idx, row in this_int_package.iterrows()]
            # flatten
            this_int_list = list(itertools.chain.from_iterable(this_int_list))
            full_int_list.append(this_int_list)

        # run this intervention package on each burnin simulation
        full_sim_list = get_combos_and_flatten([from_burnin_list, full_int_list])
        builder = ModBuilder.from_list(full_sim_list)


    ## Submit simulations ---------------------------------------------------------------------------------------
    print("Submitting")

    run_sim_args = {"config_builder": cb,
                        "exp_name": experiment_name,
                        "exp_builder": builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

