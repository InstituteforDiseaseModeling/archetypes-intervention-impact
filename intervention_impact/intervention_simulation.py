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
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report

# setup
location = 'HPC'
SetupParser.default_block = location
years = 3
max_larval_capacity = 4e8
interventions = []

sites = pd.read_csv("site_details.csv")

# Serialization
serialize = False  # If true, save serialized files
pull_from_serialization =  True # requires experiment id

for site_idx, site_name in enumerate(sites['name']):

    print("RUNNING SWEEP FOR " + site_name)

    this_site = sites.query('name==@site_name')
    exp_name = '{name}_new_site_interventions'.format(name=site_name)  # change this to something unique every time

    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                        Simulation_Duration=int(365*years),
                                        Config_Name=exp_name,
                                        Demographics_Filenames= ['demog/demog_{site}.json'.format(site=site_name)],
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

                                        # climate files
                                        Air_Temperature_Filename= os.path.join('climate', site_name,
                                                                               'air_temperature_daily.bin'),
                                        Land_Temperature_Filename= os.path.join('climate', site_name,
                                                                               'air_temperature_daily.bin'),
                                        Rainfall_Filename=os.path.join('climate', site_name,
                                                                               'rainfall_daily.bin'),
                                        Relative_Humidity_Filename=os.path.join('climate', site_name,
                                                                               'relative_humidity_daily.bin'),

                                        )

    if serialize:
        cb.update_params({'Serialization_Time_Steps': [365*years]})

    # reporting
    add_summary_report(cb)

    ## larval habitat

    all_modifications = {'arabiensis': {'habitat_split': {
                                                        'CONSTANT': 0.1,
                                                        'TEMPORARY_RAINFALL': 0.9
                                                          },
                                        'Anthropophily': 0.75
                                        },
                         'darlingi': {'habitat_split': {
                                                      "WATER_VEGETATION": 0.9,
                                                      "CONSTANT": 0.09,
                                                      "BRACKISH_SWAMP": 0.01
                                                     },
                                      'Anthropophily': 0.5
                                      },
                         'funestus': {'habitat_split': {
                                                        'WATER_VEGETATION': 1
                                                        },
                                      'Anthropophily': 0.65,
                                      'Indoor_Feeding_Fraction': 0.85

                                      },
                         'gambiae': {'habitat_split': {
                                                        'CONSTANT': 0.1,
                                                        'TEMPORARY_RAINFALL': 0.9
                                                       },
                                     'Anthropophily': 0.85,
                                     'Indoor_Feeding_Fraction': 0.85
                                     },
                         'maculatus': {'habitat_split': {
                                                        "TEMPORARY_RAINFALL": 0.8,
                                                        "WATER_VEGETATION": 0.1,
                                                        "CONSTANT": 0.1
                                                        },
                                       'Anthropophily': 0.5
                                       },
                         'minimus': {'habitat_split': {
                                                       "WATER_VEGETATION":  0.9,
                                                       "CONSTANT": 0.1
                                                      }
                                     }
                         }

    idx, vectors = this_site[[x for x in this_site.columns if x in all_modifications.keys()]].transpose().to_dict().popitem()
    vectors = {name: prop for name, prop in vectors.items() if prop>0}

    set_params_by_species(cb.params, [name for name in vectors.keys()])

    for species_name, species_prop in vectors.items():
        print("setting params for " + species_name)
        set_species_param(cb, species_name, 'Adult_Life_Expectancy', 20)

        species_mods = all_modifications[species_name]

        for param, val in species_mods.items():
            if param=='habitat_split':
                print('habitat')
                new_vals = {hab: hab_prop * species_prop * max_larval_capacity for hab, hab_prop in val.items()}
                set_species_param(cb, species_name, 'Larval_Habitat_Types', new_vals)
            else:
                set_species_param(cb, species_name, param, val)

    # itns
    def add_annual_itns(cb, year_count=1, n_rounds=1, coverage=0.8, discard_halflife=270, start_day=0):

        # per-round coverage: 1 minus the nth root of *not* getting a net in any one round
        per_round_coverage = 1 - (1 - coverage) ** (1 / n_rounds)

        for year in range(year_count):
            for round in range(n_rounds):

                add_ITN_age_season(cb,
                                   coverage_all = per_round_coverage,
                                   discard = {'halflife': discard_halflife},
                                   start=(365*year)+(30*round)+start_day)

        return {'ITN_Coverage': coverage, 'ITN_Halflife': discard_halflife, 'ITN_Per_Round_Coverage': per_round_coverage,
                'ITN_start': start_day, 'ITN_Rounds': n_rounds}

    if "itn" in interventions:
        add_annual_itns(cb, year_count=years, n_rounds=1, discard_halflife=180)


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
                                     'rate': 0.15}],
                           drug=['Artemether', 'Lumefantrine'],
                           dosing='FullTreatmentNewDetectionTech',
                           nodes={"class": "NodeSetAll"},
                           repetitions=1,
                           tsteps_btwn_repetitions=365,
                           broadcast_event_name='Received_Treatment')

        return {'ACT_Coverage': coverage}



    if "act" in interventions:
        add_healthseeking_by_coverage(cb, 0.8)


    if pull_from_serialization:
        COMPS_login("https://comps.idmod.org")
        expt = ExperimentDataStore.get_most_recent_experiment(this_site['burnin_id'].iloc[0])

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df['outpath'] = pd.Series([sim.get_path() for sim in expt.simulations])

        # temp to reduce dimensionality
        # df = df.query('x_Temporary_Larval_Habitat<100')
        # df = df.query('Run_Number==0')

        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', os.path.join(df['outpath'][x], 'output')),
            ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Filenames',
                  [name for name in os.listdir(os.path.join(df['outpath'][x], 'output')) if 'state' in name]  ),
            ModFn(DTKConfigBuilder.set_param, 'Run_Number', df['Run_Number'][x]),
            ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', df['x_Temporary_Larval_Habitat'][x]),
            ModFn(add_annual_itns, year_count=years, n_rounds=1, coverage=y/100, discard_halflife=180),
            ModFn(add_healthseeking_by_coverage,coverage=z/100),
            ModFn(add_irs_group, coverage=q/100, decay=180, start_days = [365*start for start in range(years)])
                                         ]
            for x in df.index
            for y in [0, 80]
            for z in [0, 80]
            for q in [0, 80]
        ])
    else:
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', 10**x),
            ModFn(DTKConfigBuilder.set_param, 'Run_Number', y),
            ]
            for x in np.concatenate((np.arange(0, 2.25, 0.05), np.arange(2.25, 4.25, 0.25)))
            for y in range(5)
        ])

    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    if __name__ == "__main__":
        if site_idx==0:
            SetupParser.init()
        exp_manager = ExperimentManagerFactory.init()
        exp_manager.run_simulations(**run_sim_args)