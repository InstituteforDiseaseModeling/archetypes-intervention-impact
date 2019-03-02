
import numpy as np
import itertools

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report  # add_event_counter_report
# from dtk.utils.reports.VectorReport import add_vector_stats_report
from simtools.Utilities.Experiments import retrieve_experiment

from generate_input_files import generate_input_files
from malaria.interventions.malaria_vaccine import add_vaccine
from sweep_functions import *

# variables
run_type = "intervention"  # set to "burnin" or "intervention"

# below: burnin and asset exp ids for megatrends
burnin_id = "96e9c858-a8ce-e811-a2bd-c4346bcb1555"
asset_exp_id = "96e9c858-a8ce-e811-a2bd-c4346bcb1555"

sim_root_name = "Sweep_Novel"
baseline_interventions = ["itn", "irs", "al_cm"]
baseline_intervention_coverages = [80]
sweep_interventions = ["mAb", "tbv", "pev", "atsb", "larvicides", "ivermectin"]
sweep_intervention_coverages = [0, 40, 80]
sweep_intervention_class = "list"  # can be "combo" (combinatoric) or "list"
vaccine_durations = [182, 365]
ivermectin_durations = [7, 14, 30]
new_inputs = False

# Serialization
print("setting up")
if run_type == "burnin":
    years = 15
    sweep_name = "MAP_" + sim_root_name + "_Burnin"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    sweep_name = "MAP_" + sim_root_name + "_Intervention"
    serialize = False
    pull_from_serialization = True
else:
    raise ValueError("Unknown run type " + run_type)

# setup
location = "HPC"
SetupParser.default_block = location


cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                    Simulation_Duration=int(365*years),
                                    Config_Name=sweep_name,
                                    Birth_Rate_Dependence="FIXED_BIRTH_RATE",
                                    Age_Initialization_Distribution_Type="DISTRIBUTION_COMPLEX",
                                    Num_Cores=1,

                                    # interventions
                                    # todo: do I need listed events?
                                    Listed_Events=["Bednet_Discarded", "Bednet_Got_New_One",
                                                   "Bednet_Using", "Received_Vaccine"],
                                    Enable_Default_Reporting=0,
                                    Enable_Demographics_Risk=1,
                                    Enable_Vector_Species_Report=0,

                                    # ento from prashanth
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

cb.update_params({"Disable_IP_Whitelist": 1,
                  "Enable_Property_Output": 0})

if serialize:
    cb.update_params({"Serialization_Time_Steps": [365*years]})


if __name__ == "__main__":

    SetupParser.init()

    # collect site-specific data to pass to builder functions
    COMPS_login("https://comps.idmod.org")
    sites = pd.read_csv("site_details.csv")

    print("finding collection ids and vector details")
    site_input_dir = os.path.join("sites", "all")

    with open("species_details.json") as f:
        species_details = json.loads(f.read())

    if asset_exp_id:
        print("retrieving asset experiment")
        asset_expt = retrieve_experiment(asset_exp_id)
        template_asset = asset_expt.simulations[0].tags
        cb.set_exe_collection(template_asset["exe_collection_id"])
        cb.set_dll_collection(template_asset["dll_collection_id"])
        cb.set_input_collection(template_asset["input_collection_id"])

    if new_inputs:
        print("generating input files")
        generate_input_files(site_input_dir, pop=2000, overwrite=True)

    # Find vector proportions for each vector in our site
    site_vectors = pd.read_csv(os.path.join(site_input_dir, "vector_proportions.csv"))
    simulation_setup(cb, species_details, site_vectors, site_dir=site_input_dir)

    # reporting
    for idx, row in site_vectors.iterrows():
        add_summary_report(cb,
                           age_bins=list(range(10, 130, 10)),
                           nodes={
                               "class": "NodeSetNodeList",
                               "Node_List": [int(row["node_id"])]
                           },
                           description=row["name"])
    # add_event_counter_report(cb, ["Bednet_Using", "Received_Vaccine"])
    # add_vector_stats_report(cb)

    if pull_from_serialization:
        print("building from pickup")

        # serialization
        print("retrieving burnin")
        expt = retrieve_experiment(burnin_id)

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        # temp for testing
        # df = df.iloc[0:3]

        # wrap each modfn in a tuple so it can be combined with the int lists using itertools
        from_burnin_list = [
            [ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames":
                    [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name],
                "Run_Number": df["Run_Number"][x],
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]})]
            for x in df.index]

        intervention_dict = {
            cov: {
                "itn": [ModFn(add_annual_itns,
                              year_count=years,
                              n_rounds=1,
                              coverage=cov / 100,
                              discard_halflife=180,
                              start_day=5,
                              IP=[{"NetUsage": "LovesNets"}]
                              )],
                "irs": [ModFn(add_irs_group,
                              coverage=cov / 100,
                              decay=180,
                              start_days=[365 * start for start in range(years)]
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
                                 coverage=cov/100
                                 )],
                "mAb": [ModFn(add_vaccine,
                              vaccine_type="PEV",
                              coverage=cov/100,
                              vaccine_params={
                                  "Waning_Config": {
                                      "class": "WaningEffectBox",
                                      "Box_Duration": 90
                                  }
                              },
                              target_group={"agemin": 15, "agemax": 49},
                              start_days=[365 * start for start in range(years)]
                              )],
                "pev": [ModFn(add_vaccine,
                              vaccine_type="PEV",
                              coverage=cov/100,
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
                              coverage=cov/100,
                              vaccine_params={
                                  "Reduced_Transmit": 0.75,
                                  "Waning_Config": {
                                      "class": "WaningEffectExponential",
                                      "Decay_Time_Constant": vaccine_hl / math.log(2)
                                  }
                              },
                              target_group={"agemin": 15, "agemax": 49},
                              start_days=[365 * start for start in range(years)]
                              )
                        for vaccine_hl in vaccine_durations
                        ],
                "atsb": [ModFn(add_atsb,
                               coverage=cov / 100,
                               start_days=[365 * start for start in range(years)]
                               )],
                "ors": [ModFn(add_ors,
                              coverage=cov / 100,
                              start_days=[365 * start for start in range(years)]
                              )],
                "larvicides": [ModFn(add_larvicide_wrapper,
                                     coverage=cov / 100,
                                     start_days=[365 * start for start in range(years)]
                                     )],
                "ivermectin": [ModFn(add_ivermectin_wrapper,
                                     coverage=cov / 100,
                                     drug_duration=duration,
                                     start_days=[365 * start for start in range(years)]
                                     )
                               for duration in ivermectin_durations
                               ]

            }
            for cov in list(set(baseline_intervention_coverages + sweep_intervention_coverages))
            }

        def get_combos_and_flatten(old_list):
            combos = list(itertools.product(*old_list))
            return [list(itertools.chain.from_iterable(sublist)) for sublist in combos]

        # find list of "base intervention" modfns
        baseline_int_list = [[intervention_dict[cov][intname] for cov in baseline_intervention_coverages]
                             for intname in baseline_interventions]
        baseline_int_list = get_combos_and_flatten(baseline_int_list)

        # find list of "sweep" modfns:
        if sweep_intervention_class == "list":
            sweep_int_list = [intervention_dict[cov][intname] for intname in sweep_interventions
                              for cov in sweep_intervention_coverages]
            # flatten list
            sweep_int_list = list(itertools.chain.from_iterable(sweep_int_list))

            # put everything back in a one-element list so chain.from_iterable will work
            sweep_int_list = [[item] for item in sweep_int_list]

        else:
            sweep_int_list = [[intervention_dict[cov][intname] for cov in sweep_intervention_coverages]
                              for intname in sweep_interventions]
            sweep_int_list = get_combos_and_flatten(sweep_int_list)

        full_sim_list = [from_burnin_list, baseline_int_list, sweep_int_list]
        full_sim_list = get_combos_and_flatten(full_sim_list)

        builder = ModBuilder.from_list(full_sim_list)

    else:
        print("building burnin")
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Run_Number": run_num,
                "x_Temporary_Larval_Habitat": 10 ** hab_exp}),
        ]
            for run_num in range(10)
            for hab_exp in np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))
            # for run_num in [0]
            # for hab_exp in [0, 1, 2]
        ])

    run_sim_args = {"config_builder": cb,
                    "exp_name": sweep_name,
                    "exp_builder": builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)
