import os
import pandas as pd
import sys
import json
import pdb
import re
import shutil

from dtk.vector.species import set_species_param, set_larval_habitat
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
species = 'gambiae'

main_dir = os.path.join(os.path.expanduser("~"), "Dropbox (IDM)/Malaria Team Folder/projects/Mozambique/")
hf_dir = os.path.join(main_dir, "gridded_simulation_input")
surface_dir = os.path.join(main_dir, "incidence_calibration/surface_calib/")

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                    Simulation_Duration=int(365*years),
                                    Config_Name=exp_name,
                                    Demographics_Filenames=['single_node_demographics.json'],
                                    Birth_Rate_Dependence='FIXED_BIRTH_RATE',
                                    Vector_Species_Names=[species],
                                    Num_Cores=1,

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

set_species_param(cb, species, "Adult_Life_Expectancy", 20)

ento_trend =   [0.0429944166751962,
               0.145106159922212,
               0.220520011001099,
               0.318489404300663,
               0.0617610600835594,
               0.0462380862878181,
               0.0367590381502996,
               0.02474944109524821,
               0.0300445801767523,
               0.021859890543704,
               0.0261404367939001,
               0.0253992634551118]
hab = {species : {
        # 'CONSTANT': 2e6,
        "LINEAR_SPLINE": {
                           "Capacity_Distribution_Per_Year": {
                               "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083, 182.5, 212.917, 243.333, 273.75,
                                         304.167, 334.583],
                               "Values": ento_trend
                           },
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

