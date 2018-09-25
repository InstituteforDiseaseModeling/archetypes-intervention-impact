import pandas as pd
import json

import os

from dtk.vector.species import set_params_by_species, set_species_param
from dtk.interventions.health_seeking import add_health_seeking
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.property_change import change_individual_property

def assign_net_ip(cb, hates_net_prop):
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop),
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop,
          trigger_condition_list=["Births"])
    return {"Hates_Nets": hates_net_prop}

def simulation_setup(cb, species_details, max_larval_capacity=4e8):

    site_dir = os.path.join("sites", "all")

    # directories
    cb.update_params({
                    "Demographics_Filenames": [os.path.join(site_dir, "demographics.json"),
                                               os.path.join(site_dir, "demographics_net_overlay.json")],
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

    # Find vector proportions for each vector TODO: add node-specific and habitat-specific multiplier to demog overlay
    set_params_by_species(cb.params, [name for name in species_details.keys()])

    for species_name, species_modifications in species_details.items():
        set_species_param(cb, species_name, "Adult_Life_Expectancy", 20)
        set_species_param(cb, species_name, "Vector_Sugar_Feeding_Frequency", "VECTOR_SUGAR_FEEDING_NONE")

        for param, val in species_modifications.items():
            if param == "habitat_split":
                new_vals = {hab: hab_prop * max_larval_capacity for hab, hab_prop in val.items()}
                set_species_param(cb, species_name, "Larval_Habitat_Types", new_vals)
            else:
                set_species_param(cb, species_name, param, val)



# itns
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
            "ITN_start": start_day, "ITN_Rounds": n_rounds, "ITN_Distributions": year_count}


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
        add_IRS(cb, start, [{"min": 0, "max": 200, "coverage": coverage}],
                waning=waning)

    return {"IRS_Halflife": decay, "IRS_Start": start_days[0], "IRS_Coverage": coverage}


# act
def add_healthseeking_by_coverage(cb, coverage=1.0, rate=0.15):
    add_health_seeking(cb,
                       targets=[{"trigger": "NewClinicalCase",
                                 "coverage": coverage,
                                 "agemin": 0,
                                 "agemax": 100,
                                 "seek": 1.0,
                                 "rate": rate}],
                       drug=["Artemether", "Lumefantrine"],
                       dosing="FullTreatmentNewDetectionTech",
                       nodes={"class": "NodeSetAll"},
                       repetitions=1,
                       tsteps_btwn_repetitions=365,
                       broadcast_event_name="Received_Treatment")

    return {"ACT_Coverage": coverage, "ACT_Daily_Prob": rate}
