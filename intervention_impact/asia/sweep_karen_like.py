import os
import pandas as pd
import sys
import json
import pdb
import re
import shutil

from dtk.vector.species import set_species_param, set_larval_habitat, set_params_by_species
from simtools.SetupParser import SetupParser
from dtk.interventions.health_seeking import add_health_seeking
from dtk.generic.climate import set_climate_constant
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_event_counter_report

# setup
location = 'HPC'
SetupParser.default_block = location
exp_name = 'Intervention_Burnin_Test'  # change this to something unique every time
years = 50
species = 'minimus'

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

builder = ModBuilder.from_list([[
    ModFn(DTKConfigBuilder.set_param, 'Run_Number', i),
    ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', x
          )
    ]
    for i, x in enumerate([5, 10, 50, 100, 500, 1000])
])

run_sim_args = {'config_builder': cb,
                'exp_name': exp_name,
                'exp_builder': builder}

