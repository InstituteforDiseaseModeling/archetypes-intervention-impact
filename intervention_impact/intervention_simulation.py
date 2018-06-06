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
from dtk.interventions.irs import add_IRS
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
archetype = "moine"
exp_name = 'Moine_ITN_CM_80_Climate'  # change this to something unique every time
years = 10
interventions = ['itn']
use_climate = True

# Serialization
serialize = False  # If true, save serialized files
pull_from_serialization =  True # requires experiment id

archetypes = {'karen': {
                        'demog': 'demog/demog_karen.json',
                        'species': [{'name': 'minimus',
                                    'seasonality': {
                                      "Times":  [0, 1, 244, 274, 363],
                                      "Values": [0.2, 0.2, 0.7, 3, 3]
                                      }
                                    }],
                        'burnin_id': "e88827cf-a365-e811-a2c0-c4346bcb7275"
             },
             'moine': {
                        'demog': 'demog/demog_moine.json',
                        'species': [{'name':'gambiae',
                                     'seasonality': {
                                         "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083, 182.5, 212.917,
                                                   243.333, 273.75,
                                                   304.167, 334.583],
                                         "Values": [0.0429944166751962,
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
                                     }
                        }],
                        'burnin_id': "476808e8-6369-e811-a2c0-c4346bcb7275"

             }

}

arch_vals = archetypes[archetype]

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                    Simulation_Duration=int(365*years),
                                    Config_Name=exp_name,
                                    Demographics_Filenames=[arch_vals['demog']],
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

if use_climate:
    climate_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'data',
                               'Mozambique', 'Magude', 'Mozambique climate')
    cb.update_params( # Climate file paths
        {'Air_Temperature_Filename': os.path.join('climate',
                                                 'Mozambique_Moine_2.5arcmin_air_temperature_daily.bin'),
        'Land_Temperature_Filename': os.path.join('climate',
                                                  'Mozambique_Moine_2.5arcmin_air_temperature_daily.bin'),
        'Rainfall_Filename': os.path.join('climate',
                                          'Mozambique_Moine_2.5arcmin_rainfall_daily.bin'),
        'Relative_Humidity_Filename': os.path.join('climate',
                                                   'Mozambique_Moine_2.5arcmin_relative_humidity_daily.bin')})

    total_capacity = 4e8

    set_params_by_species(cb.params, ['funestus', 'gambiae'])
    set_larval_habitat(cb, {'funestus': {'WATER_VEGETATION': total_capacity*0.94},
                            'gambiae': {'CONSTANT': total_capacity*0.06*0.9,
                                        'TEMPORARY_RAINFALL': total_capacity*0.06*0.1}})

else:
    set_climate_constant(cb)
    set_params_by_species(cb.params, [species['name'] for species in arch_vals['species']])

    for species in arch_vals['species']:

        print('setting params for species ' + species['name'])

        set_species_param(cb, species['name'], "Adult_Life_Expectancy", 20)

        hab = {species['name'] : {
                "LINEAR_SPLINE": {
                                   "Capacity_Distribution_Over_Time": species['seasonality'],
                                   "Capacity_Distribution_Number_Of_Years": 1,
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

# itns
def add_annual_itns(cb, year_count=1, coverage=0.8, discard_halflife=270):
    for year in range(year_count):
        add_ITN_age_season(cb,
                           coverage_all = coverage,
                           discard = {'halflife': discard_halflife},
                           start=365*year)

    return {'ITN_Coverage': coverage, 'ITN_Halflife': discard_halflife}

if "itn" in interventions:
    add_annual_itns(cb, year_count=years)


# irs
def add_irs_group(cb, coverage=1.0, start_days=[0], decay=270):

    waning = {
                "Killing_Config": {
                    "Initial_Effect": 0.6,
                    "Decay_Time_Constant": decay,
                    "class": "WaningEffectExponential"
                },
                "Blocking_Config": {
                    "Initial_Effect": 0.0,
                    "Decay_Time_Constant": 730,
                    "class": "WaningEffectExponential"
                }}

    for start in start_days:
        add_IRS(cb, start, [{'min': 0, 'max': 200, 'coverage': coverage}],
                waning=waning)

    return {'IRS_Halflife': decay, 'IRS_Start': start_days[0], 'IRS_Coverage': coverage}

if "irs" in interventions:
    add_irs_group(cb, coverage=0.8, decay=180) # simulate actellic irs

# act
def add_healthseeking_by_coverage(cb, coverage=1.0):
    add_health_seeking(cb,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': coverage, 'agemin': 0, 'agemax': 100, 'seek': 1.0,
                                 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'],
                       dosing='FullTreatmentNewDetectionTech',
                       nodes={"class": "NodeSetAll"},
                       repetitions=1,
                       tsteps_btwn_repetitions=365,
                       broadcast_event_name='Received_Treatment')

    return {'Healthseek_Coverage': coverage}



if "act" in interventions:
    add_healthseeking_by_coverage(cb, 0.8)


if pull_from_serialization:
    COMPS_login("https://comps.idmod.org")
    expt = ExperimentDataStore.get_most_recent_experiment(arch_vals['burnin_id'])

    df = pd.DataFrame([x.tags for x in expt.simulations])
    df['outpath'] = pd.Series([sim.get_path() for sim in expt.simulations])

    # temp to reduce dimensionality
    # df = df.query('x_Temporary_Larval_Habitat==100')

    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', os.path.join(df['outpath'][x], 'output')),
        ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Filenames',
              [name for name in os.listdir(os.path.join(df['outpath'][x], 'output')) if 'state' in name]  ),
        ModFn(DTKConfigBuilder.set_param, 'Run_Number', df['Run_Number'][x]),
        ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', df['x_Temporary_Larval_Habitat'][x]),
        # ModFn(add_annual_itns, year_count=years, coverage=z/100),
        # ModFn(add_healthseeking_by_coverage,coverage=zz/100),
        # ModFn(add_irs_group, coverage=z/100, decay=180)
                                    ]
        for x in df.index # for z in range(0, 100, 20)  # for zz in range(0, 100, 20)
    ])
else:
    builder = ModBuilder.from_list([[
        ModFn(DTKConfigBuilder.set_param, 'Run_Number', y),
        ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', 10**x),
        # ModFn(add_ITN_age_season, start=365*2, coverage_all=z/100)
        ]
        for x in np.concatenate((np.arange(0, 2.25, 0.05), np.arange(2.25, 4.25, 0.25))) for y in range(5) # for z in [0, 50, 80]
        # for x in [0.5, 0, 1] # for y in range(2)
    ])

run_sim_args = {'config_builder': cb,
                'exp_name': exp_name,
                'exp_builder': builder}

