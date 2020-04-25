###############################################################################################################
## intervention_functions.py
## Amelia Bertozzi-Villa
## October 2019
##
## Library of functions to set up sweeps of different interventions.
##############################################################################################################

## Importing  ---------------------------------------------------------------------------------------

import pandas as pd
import json
import math

import os
import pdb

from simtools.Utilities.Experiments import retrieve_experiment
from simtools.ModBuilder import ModBuilder, ModFn

from dtk.vector.species import set_params_by_species, set_species_param
from dtk.interventions.habitat_scale import scale_larval_habitats
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.property_change import change_individual_property
from dtk.interventions.novel_vector_control import add_ATSB, add_ors_node, add_larvicides
from dtk.interventions.ivermectin import add_ivermectin
from malaria.interventions.malaria_vaccine import add_vaccine
from malaria.interventions.adherent_drug import configure_adherent_drug

from malaria.reports.MalariaReport import add_summary_report
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign, add_diagnostic_survey
from malaria.interventions.health_seeking import add_health_seeking

## Interventions ---------------------------------------------------------------------------------------

def generate_intervention_tuples(coverages, start_days, years,
                                 vaccine_durations=None, ivermectin_durations=None, smc_max_ages=None):

    # defaults for durations until we find a better way to pass these
    if ivermectin_durations is None:
        ivermectin_durations = [14]
    if vaccine_durations is None:
        vaccine_durations = [365]
    if smc_max_ages is None:
        smc_max_ages = [5]
    intervention_dict = {
        start_day: {
            cov: {
                "itn": [ModFn(add_annual_itns,
                              year_count=years,
                              n_rounds=1,
                              coverage=cov / 100,
                              discard_halflife=180,
                              start_day=start_day,
                              IP=[{"NetUsage": "LovesNets"}]
                              )],
                "irs": [ModFn(add_irs_group,
                              coverage=cov / 100,
                              decay=180,
                              start_days=[(365 * year_idx) + start_day for year_idx in range(years)]
                              )],
                "al_cm": [ModFn(add_healthseeking_by_coverage,
                                coverage=cov / 100,
                                rate=0.15
                                )],
                "dp_cm": [ModFn(add_healthseeking_by_coverage,
                                coverage=cov / 100,
                                rate=0.15,
                                drugname="DP"
                                )],
                "dp_mda": [ModFn(add_mda,
                                 coverage=cov / 100
                                 )],
                "mAb": [ModFn(add_vaccine,
                              vaccine_type="PEV",
                              coverage=cov / 100,
                              vaccine_params={
                                  "Waning_Config": {
                                      "class": "WaningEffectBox",
                                      "Box_Duration": 90
                                  }
                              },
                              target_group={"agemin": 15, "agemax": 49},
                              start_days=[(365 * year_idx) + start_day for year_idx in range(years)]
                              )],
                "pev": [ModFn(add_vaccine,
                              vaccine_type="PEV",
                              coverage=cov / 100,
                              vaccine_params={
                                  "Reduced_Acquire": 0.75,
                                  "Waning_Config": {
                                      "class": "WaningEffectExponential",
                                      "Decay_Time_Constant": vaccine_hl / math.log(2)
                                  }
                              },
                              trigger_condition_list=["Births"],
                              start_days=[0],
                              triggered_delay=182)
                        for vaccine_hl in vaccine_durations
                        ],
                "tbv": [ModFn(add_vaccine,
                              vaccine_type="TBV",
                              coverage=cov / 100,
                              vaccine_params={
                                  "Reduced_Transmit": 0.75,
                                  "Waning_Config": {
                                      "class": "WaningEffectExponential",
                                      "Decay_Time_Constant": vaccine_hl / math.log(2)
                                  }
                              },
                              target_group={"agemin": 15, "agemax": 49},
                              start_days=[(365 * year_idx) + start_day for year_idx in range(years)]
                              )
                        for vaccine_hl in vaccine_durations
                        ],
                "atsb": [ModFn(add_atsb,
                               coverage=cov / 100,
                               # renew atsb every 6 months
                               start_days=[(365 * year_idx) + start_day + midyear_start
                                           for year_idx in range(years) for midyear_start in [0, 182]],
                               initial_effect= 1
                               )
                         ],
                "ors": [ModFn(add_ors,
                              coverage=cov / 100,
                              start_days=[(365 * year_idx) + start_day for year_idx in range(years)]
                              )],
                "larvicides": [ModFn(add_larvicide_wrapper,
                                     coverage=cov / 100,
                                     start_days=[(365 * year_idx) + start_day for year_idx in range(years)]
                                     )],
                "ivermectin": [ModFn(add_ivermectin_wrapper,
                                     coverage=cov / 100,
                                     drug_duration=duration,
                                     start_days=[(365 * year_idx) + start_day for year_idx in range(years)],
                                     monthly_rounds=3
                                     )
                               for duration in ivermectin_durations
                               ],
                "smc": {max_age: [ModFn(add_smc,
                              coverage=cov/100,
                              start_day=start_day,
                              agemax=max_age
                              )]
                        for max_age in smc_max_ages
                        }

            }
            for cov in coverages
        }
        for start_day in start_days
    }

    return(intervention_dict)


def add_annual_itns(cb, year_count=1, n_rounds=1, coverage=0.8, discard_halflife=270, start_day=0, IP=[]):

    # per-round coverage: 1 minus the nth root of *not* getting a net in any one round
    per_round_coverage = 1 - (1 - coverage) ** (1 / n_rounds)

    for year in range(year_count):
        for this_round in range(n_rounds):

            add_ITN_age_season(cb,
                               coverage_all=per_round_coverage,
                               discard={"halflife": discard_halflife},
                               start=(365 * year) + (30 * this_round) + start_day,
                               ind_property_restrictions=IP)

    return {"ITN_Coverage": coverage, "ITN_Halflife": discard_halflife, "ITN_Per_Round_Coverage": per_round_coverage,
            "ITN_Start": start_day, "ITN_Rounds": n_rounds, "ITN_Distributions": year_count}


def add_annual_itns_w_irs(cb, year_count=1, n_rounds=1, coverage=0.8, discard_halflife=270, start_day=0, IP=[]):

    # per-round coverage: 1 minus the nth root of *not* getting a net in any one round
    per_round_coverage = 1 - (1 - coverage) ** (1 / n_rounds)

    for year in range(year_count):
        for this_round in range(n_rounds):

            add_ITN_age_season(cb,
                               coverage_all=per_round_coverage,
                               discard={"halflife": discard_halflife},
                               start=(365 * year) + (30 * this_round) + start_day,
                               ind_property_restrictions=IP)

    add_irs_group(cb, decay=180, trigger_condition_list=['Bednet_Got_New_One'])
    return {"ITN_IRS_Coverage": coverage, "ITN_Halflife": discard_halflife, "ITN_Per_Round_Coverage": per_round_coverage,
            "ITN_start": start_day, "ITN_Rounds": n_rounds, "ITN_Distributions": year_count}


# irs
def add_irs_group(cb, coverage=1.0, start_days=[0], decay=270, trigger_condition_list=[]):

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
        add_IRS(cb, start, [{"min": 0, "max": 200, "coverage": coverage}],
                waning=waning, trigger_condition_list=trigger_condition_list)

    return {"IRS_Halflife": decay, "IRS_Start": start_days[0], "IRS_Coverage": coverage}


# act
def add_healthseeking_by_coverage(cb, coverage=1.0, rate=0.15, drugname="AL"):
    drugs = {"AL": ["Artemether", "Lumefantrine"],
             "DP": ["DHA", "Piperaquine"]}

    add_health_seeking(cb,
                       targets=[{"trigger": "NewClinicalCase",
                                 "coverage": coverage,
                                 "agemin": 0,
                                 "agemax": 100,
                                 "seek": 1.0,
                                 "rate": rate}],
                       drug=drugs[drugname],
                       dosing="FullTreatmentNewDetectionTech",
                       nodes=None,
                       repetitions=1,
                       tsteps_btwn_repetitions=365,
                       broadcast_event_name="Received_Treatment")

    return {"CM_Coverage": coverage, "CM_Daily_Prob": rate, "CM_Drug": drugname}

# mda
def add_mda(cb, coverage=0.8, drugname="DP", start_days=[0], reps=3):

    add_drug_campaign(cb, 'MDA', drugname,
                      start_days=start_days,
                      coverage=coverage,
                      # nodes={"class": "NodeSetAll"},
                      repetitions=reps,
                      tsteps_btwn_repetitions=30)

    return {"MDA_Drug": drugname, "MDA_Repetitions": reps, "MDA_Coverage": coverage}


def add_atsb(cb, coverage=0.8, start_days=[0], duration=180, initial_effect=0.0377):

    for start_day in start_days:

        # default killing cfg
        orig_killing_cfg = {"Killing_Config": {
            "class": "WaningEffectBoxExponential",
            "Box_Duration": 180,
            "Decay_Time_Constant": 130,
            "Initial_Effect": initial_effect
        }}

        add_ATSB(cb, coverage=coverage,
                 start_day=start_day,
                 duration=duration,
                 kill_cfg=orig_killing_cfg
                 )

    return {"ATSB_Coverage": coverage, "ATSB_Start": start_days[0], "ATSB_Duration_Mean": duration,
            "ATSB_Initial_Effect": initial_effect}


def add_ors(cb, coverage, start_days=[0], kill_duration=100):

    for start_day in start_days:
        add_ors_node(cb, start_day=start_day,
                     coverage=coverage,
                     killing_decay=kill_duration)

    return {"ORS_Coverage": coverage, "ORS_Start": start_days[0], "ORS_Halflife": kill_duration}


def add_larvicide_wrapper(cb, coverage, start_days=[0], kill_duration=100):

    for start_day in start_days:
        add_larvicides(cb, start_day=start_day,
                     coverage=coverage,
                     killing_duration=kill_duration)

    return {"Larvicide_Coverage": coverage, "Larvicide_Start": start_days[0], "Larvicide_Kill_Halflife": kill_duration}

def add_ivermectin_wrapper(cb, coverage, start_days=[0], drug_duration=7, monthly_rounds=3):
    start_days = [(30 * month_idx) + start for start in start_days for month_idx in range(monthly_rounds)]

    add_ivermectin(cb, start_days=start_days,
                   coverage=coverage,
                   box_duration=drug_duration)

    return {"Ivermectin_Coverage": coverage, "Ivermectin_Start": start_days[0], "Ivermectin_Duration": drug_duration,
            "Ivermectin_Monthly_Rounds": monthly_rounds}

def smc_adherent_configuration(cb, adherence):
    smc_adherent_config = configure_adherent_drug(cb, doses = [["Sulfadoxine", "Pyrimethamine",'Amodiaquine'],
                                                         ['Amodiaquine'],
                                                         ['Amodiaquine']],
                                                  dose_interval=1,
                                                  non_adherence_options=['Stop'],
                                                  non_adherence_distribution=[1],
                                                  adherence_config={
                                                        "class": "WaningEffectMapCount",
                                                        "Initial_Effect": 1,
                                                        "Durability_Map": {
                                                            "Times": [
                                                                1.0,
                                                                2.0,
                                                                3.0
                                                            ],
                                                            "Values": [
                                                                1, # for day 1
                                                                adherence, # day 2
                                                                adherence # day 3
                                                            ]
                                                        }
                                                    }
    )
    return smc_adherent_config

def add_smc(cb, coverage, start_day=0, adherence=0.6, agemax=5):

    adherent_drug_configs = smc_adherent_configuration(cb, adherence)  # adherence
    add_diagnostic_survey(cb, start_day=start_day,
                          coverage=1,
                          target={"agemin": 0.25, "agemax": agemax},
                          diagnostic_type='FEVER',
                          diagnostic_threshold=0.5,
                          negative_diagnosis_configs=[{
                              "Broadcast_Event": "No_SMC_Fever",
                              "class": "BroadcastEvent"}],
                          positive_diagnosis_configs=[{
                              "Broadcast_Event": "Has_SMC_Fever",
                              "class": "BroadcastEvent"}]
                          )
    add_drug_campaign(cb, 'SMC', 'SPA', start_days=[start_day],
                      coverage=coverage, target_group={'agemin': 0.25, 'agemax': agemax},
                      listening_duration=2,
                      trigger_condition_list=['No_SMC_Fever'],
                      adherent_drug_configs=[adherent_drug_configs])

    return {"SMC_Coverage": coverage, "SMC_Start": start_day, "SMC_Adherence": adherence,
            "SMC_Max_Age": agemax}
