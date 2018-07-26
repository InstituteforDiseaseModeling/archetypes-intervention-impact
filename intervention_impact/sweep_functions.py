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
    # change_individual_property(cb, "NetUsage", "LovesNets", coverage=1 - hates_net_prop),
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop,
          trigger_condition_list=["Births"]),
    # change_individual_property(cb, "NetUsage", "LovesNets", coverage=1 - hates_net_prop,
    #       trigger_condition_list=["Births"])

    return {"Hates_Nets": hates_net_prop}

def site_simulation_setup(cb, site_name, species_details, vectors, max_larval_capacity=4e8):

    site_dir = os.path.join("sites", site_name)

    # directories
    # todo: do I need this if I specify an asset collection?
    cb.update_params({
                    "Demographics_Filenames": ["sites/{site}/demographics_{site}.json".format(site=site_name),
                                               "sites/{site}/demographics_{site}_hatenets_0.json".format(site=site_name)],
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

    # Find vector proportions for each vector in our site
    set_params_by_species(cb.params, [name for name in vectors.keys()])

    for species_name, species_prop in vectors.items():
        set_species_param(cb, species_name, "Adult_Life_Expectancy", 20)
        set_species_param(cb, species_name, "Vector_Sugar_Feeding_Frequency", "VECTOR_SUGAR_FEEDING_NONE")

        species_modifications = species_details[species_name]

        for param, val in species_modifications.items():
            if param == "habitat_split":
                new_vals = {hab: hab_prop * species_prop * max_larval_capacity for hab, hab_prop in val.items()}
                set_species_param(cb, species_name, "Larval_Habitat_Types", new_vals)
            else:
                set_species_param(cb, species_name, param, val)

    return_params = {"{species}_proportion".format(species=species): prop for species, prop in vectors.items()}
    return_params["Site_Name"] = site_name

    return return_params


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

    return {"ACT_Coverage": coverage, "ACT_HS_Rate": rate}
