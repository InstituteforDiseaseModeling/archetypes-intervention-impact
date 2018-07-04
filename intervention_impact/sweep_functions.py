
import pandas as pd
import os
import json


from dtk.vector.species import set_params_by_species, set_species_param
from dtk.interventions.health_seeking import add_health_seeking
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season

def site_simulation_setup(cb, site_name, pull_serialized=True, max_larval_capacity=4e8):
    sites = pd.read_csv("site_details.csv")
    this_site = sites.query('name==@site_name')

    # directories
    cb.update_params({
                    "Demographics_Filenames": ["demog/demog_{site}.json".format(site=site_name)],
                    "Air_Temperature_Filename": os.path.join("climate", site_name,
                                                           "air_temperature_daily.bin"),
                    "Land_Temperature_Filename": os.path.join('climate', site_name,
                                                           "air_temperature_daily.bin"),
                    "Rainfall_Filename": os.path.join("climate", site_name,
                                                           "rainfall_daily.bin"),
                    "Relative_Humidity_Filename": os.path.join("climate", site_name,
                                                           "relative_humidity_daily.bin")
                    }
    )

    # vectors
    with open("species_details.json") as f:
        species_details = json.loads(f.read())

    # Find vector details for each vector in our site (todo: make this more interpretable-- pull from site-specific file)
    idx, vectors = this_site[[x for x in this_site.columns if x in species_details.keys()]].transpose().to_dict().popitem()
    vectors = {name: prop for name, prop in vectors.items() if prop > 0}

    set_params_by_species(cb.params, [name for name in vectors.keys()])

    for species_name, species_prop in vectors.items():
        set_species_param(cb, species_name, 'Adult_Life_Expectancy', 20)

        species_modifications = species_details[species_name]

        for param, val in species_modifications.items():
            if param=='habitat_split':
                new_vals = {hab: hab_prop * species_prop * max_larval_capacity for hab, hab_prop in val.items()}
                set_species_param(cb, species_name, 'Larval_Habitat_Types', new_vals)
            else:
                set_species_param(cb, species_name, param, val)

    return {"Site_Name": site_name, "Vectors": vectors} # todo: can comps accept a dict as a key

# itns
def add_annual_itns(cb, year_count=1, n_rounds=1, coverage=0.8, discard_halflife=270, start_day=0):

    # per-round coverage: 1 minus the nth root of *not* getting a net in any one round
    per_round_coverage = 1 - (1 - coverage) ** (1 / n_rounds)

    for year in range(year_count):
        for round in range(n_rounds):

            add_ITN_age_season(cb,
                               coverage_all = per_round_coverage,
                               discard = {'halflife': discard_halflife},
                               start=( 365 *year ) +( 30 *round ) +start_day)

    return {'ITN_Coverage': coverage, 'ITN_Halflife': discard_halflife, 'ITN_Per_Round_Coverage': per_round_coverage,
            'ITN_start': start_day, 'ITN_Rounds': n_rounds}

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

# act
def add_healthseeking_by_coverage(cb, coverage=1.0):
    add_health_seeking(cb,
                       targets=
                           [{'trigger': 'NewClinicalCase', 'coverage': coverage, 'agemin': 0, 'agemax': 100, 'seek': 1.0,
                                 'rate': 0.15}],
                       drug=['Artemether', 'Lumefantrine'],
                       dosing='FullTreatmentNewDetectionTech',
                       nodes={"class": "NodeSetAll"},
                       repetitions=1,
                       tsteps_btwn_repetitions=365,
                       broadcast_event_name='Received_Treatment')

    return {'ACT_Coverage': coverage}
