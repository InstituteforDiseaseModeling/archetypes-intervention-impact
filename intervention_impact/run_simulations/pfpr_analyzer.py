import os
import pdb

import pandas as pd
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import exps_for_suite_id


class PfPRAnalyzer(BaseAnalyzer):

    def __init__(self, dir_name, report_names=["AnnualAverage"], sweep_variables=None, working_dir="."):
        super(PfPRAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.sitenames=report_names
        self.dir_name = dir_name

    def select_simulation_data(self, data, simulation):
        colname = "initial_prev" if self.dir_name == "initial" else "final_prev"

        simdata = []

        for site_name in self.sitenames:

            try:
                channeldata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTime"]["PfPR_2to10"]
            except:
                print("file not found for sim" + simulation.id)

            tempdata = pd.DataFrame({colname: channeldata,
                                    "Site_Name": site_name})
            tempdata = tempdata[-2:-1]
            simdata.append(tempdata)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def filter(self, simulation):
        return simulation.id != "87227af9-0fcf-e811-a2bd-c4346bcb1555"

    def finalize(self, all_data):
        data_sets_per_experiment = {}

        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.exp_name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)
            d.to_csv(os.path.join(self.working_dir, self.dir_name, "{name}.csv".format(name=experiment_name)), index=False)


if __name__ == "__main__":
    SetupParser.init("HPC")
    out_dir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
                           'map_intervention_impact', 'lookup_tables')

    sites = pd.read_csv("site_details.csv")

    experiments = {# "initial": "96e9c858-a8ce-e811-a2bd-c4346bcb1555",
                   "interactions" :"636d798d-01d1-e811-a2bd-c4346bcb1555"
                   }

    for dirname, exp_id in experiments.items():

        am = AnalyzeManager(exp_list=exp_id, analyzers=[PfPRAnalyzer(working_dir=out_dir,
                                                                     dir_name=dirname,
                                                                     report_names = sites["name"].tolist(),
                                                                      sweep_variables=["Run_Number",
                                                                                       "x_Temporary_Larval_Habitat",
                                                                                       "CM_Drug",
                                                                                       "CM_Coverage",
                                                                                       "IRS_Coverage",
                                                                                       "ITN_Coverage",
                                                                                       "MDA_Drug",
                                                                                       "MDA_Repetitions",
                                                                                       "MDA_Coverage",
                                                                                       "PEV_Coverage",
                                                                                       "PEV_Waning_Config_class",
                                                                                       "PEV_Waning_Config_Decay_Time_Constant",
                                                                                       "TBV_Coverage",
                                                                                       "TBV_Waning_Config_Decay_Time_Constant",

                                                                                       ])],
                            force_analyze=True)

        print(am.experiments)
        am.analyze()

