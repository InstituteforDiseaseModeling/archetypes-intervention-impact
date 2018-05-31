import os
import pandas as pd
import numpy as np
import sys
import json
import pdb
import re
import shutil

from dtk.vector.species import set_species_param, set_larval_habitat, set_params_by_species
from simtools.SetupParser import SetupParser
from dtk.interventions.health_seeking import add_health_seeking
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.generic.climate import set_climate_constant
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report

# setup
location = 'HPC'
SetupParser.default_block = location
exp_name = 'Karen_Fine_Burnin'  # change this to something unique every time
years = 50
species = 'minimus'

# Serialization
serialize = True  # If true, save serialized files
pull_from_serialization =  False # requires serialization date and experiment id
serialization_exp_id = "1c52b2f5-4064-e811-a2c0-c4346bcb7275"

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                    Simulation_Duration=int(365*years),
                                    Config_Name=exp_name,
                                    Demographics_Filenames=['demo_karen.json'],
                                    Birth_Rate_Dependence='FIXED_BIRTH_RATE',
                                    Num_Cores=1,

                                    # interventions
                                    Valid_Intervention_States= [],  # apparently a necessary parameter
                                    Listed_Events= ['Bednet_Discarded', 'Bednet_Got_New_One', 'Bednet_Using'],

                                    ## ento from prashanth
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
    cb.update_params({'Serialization_Time_Steps': [365*years]})

# reporting
add_summary_report(cb)

## larval habitat
set_climate_constant(cb)

set_params_by_species(cb.params, ['minimus'])
set_species_param(cb, species, "Adult_Life_Expectancy", 20)

capacity_dist_per_year = {
                          "Times":  [0, 1, 244, 274, 363],
                          "Values": [0.2, 0.2, 0.7, 3, 3]
                      }
hab = {species : {
        # 'CONSTANT': 2e6,
        "LINEAR_SPLINE": {
                           "Capacity_Distribution_Per_Year": capacity_dist_per_year,
                           "Max_Larval_Capacity": 1e8
                       }
                       }
           }

set_larval_habitat(cb, hab)

cb.update_params({
        "Report_Event_Recorder": 1,
        "Report_Event_Recorder_Events": ["Bednet_Using"],
        "Report_Event_Recorder_Ignore_Events_In_List": 0
    })

if pull_from_serialization:
    COMPS_login("https://comps.idmod.org")
    expt = ExperimentDataStore.get_most_recent_experiment(serialization_exp_id)

    df = pd.DataFrame([x.tags for x in expt.simulations])
    df['outpath'] = pd.Series([sim.get_path() for sim in expt.simulations])

    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', os.path.join(df['outpath'][x], 'output')),
        ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Filenames',
              [name for name in os.listdir(os.path.join(df['outpath'][x], 'output')) if 'state' in name]  ),
        ModFn(DTKConfigBuilder.set_param, 'Run_Number', y),
        ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', df['x_Temporary_Larval_Habitat'][x]),
        ModFn(add_ITN_age_season, coverage_all=z/100)
                                    ]
        for x in df.index for y in range(10) for z in [0, 80, 100] #range(0,105, 5)
    ])
else:
    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.set_param, 'Run_Number', y),
        ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', 10**x),
        # ModFn(add_ITN_age_season, start=365*2, coverage_all=z/100)
        ]
        for x in np.arange(0, 4.25, 0.25) for y in range(10) # for z in [0, 50, 80]
        # for x in [100] for y in [0,50,80]
    ])

run_sim_args = {'config_builder': cb,
                'exp_name': exp_name,
                'exp_builder': builder}

