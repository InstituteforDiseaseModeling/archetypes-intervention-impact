import os
import pdb

import pandas as pd
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import exps_for_suite_id


class PfPRAnalyzer(BaseAnalyzer):

    def __init__(self, sweep_variables=None, working_dir="."):
        super(PfPRAnalyzer, self).__init__(working_dir=working_dir,
                                           filenames=["output/MalariaSummaryReport_AnnualAverage.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]

    def select_simulation_data(self, data, simulation):
        colname = "initial_prev" if simulation.experiment.exp_name == "initial" else "final_prev"
        age_bins = data[self.filenames[0]]["Metadata"]["Age Bins"]
        channeldata = data[self.filenames[0]]["DataByTimeAndAgeBins"]["Average Population by Age Bin"]
        simdata = pd.DataFrame(channeldata).transpose()
        simdata["age"] = age_bins
        print(simdata)


        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):
        data_sets_per_experiment = {}

        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.exp_name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)
            d.to_csv(os.path.join(self.working_dir, "{name}.csv".format(name=experiment_name)), index=False)


if __name__ == "__main__":
    SetupParser.init("HPC")
    out_dir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
                           'map_intervention_impact', 'lookup_tables')

    run_type = "exp"

    if run_type == "exp":
        experiments = ["6e91248b-dbbd-e811-a2bd-c4346bcb1555"]

        for exp_id in experiments:

            am = AnalyzeManager(exp_list=exp_id, analyzers=[PfPRAnalyzer(working_dir=out_dir,
                                                                              sweep_variables=["Site_Name",
                                                                                               "Run_Number",
                                                                                               "x_Temporary_Larval_Habitat",
                                                                                               "ACT_Coverage",
                                                                                               "IRS_Coverage",
                                                                                               "ITN_Coverage"
                                                                                               ])],
                                force_analyze=True)

            print(am.experiments)
            am.analyze()


    elif run_type == "suite":

        exps = exps_for_suite_id("537b9041-0fa3-e811-a2c0-c4346bcb7275")
        print(exps)
        am = AnalyzeManager(exp_list=exps, analyzers=[PfPRAnalyzer(working_dir=out_dir,
                                                                          sweep_variables=["Site_Name",
                                                                                           "Run_Number",
                                                                                           "x_Temporary_Larval_Habitat",
                                                                                           "ACT_Coverage",
                                                                                           "IRS_Coverage",
                                                                                           "ITN_Coverage"
                                                                                           ])])
        am.analyze()
