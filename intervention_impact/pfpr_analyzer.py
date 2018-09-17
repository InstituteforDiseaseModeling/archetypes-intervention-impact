import os

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
        channeldata = data[self.filenames[0]]["DataByTime"]["PfPR_2to10"]
        simdata = pd.DataFrame({colname: channeldata})
        simdata = simdata[-2:-1]

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
        experiments = ["547b9041-0fa3-e811-a2c0-c4346bcb7275", "557b9041-0fa3-e811-a2c0-c4346bcb7275",
                       "567b9041-0fa3-e811-a2c0-c4346bcb7275", "577b9041-0fa3-e811-a2c0-c4346bcb7275",
                       "587b9041-0fa3-e811-a2c0-c4346bcb7275", "5a7b9041-0fa3-e811-a2c0-c4346bcb7275",
                       "5b7b9041-0fa3-e811-a2c0-c4346bcb7275", "c6ddc170-19a1-e811-a2c0-c4346bcb7275"]

        for exp_id in experiments:

            am = AnalyzeManager(exp_list=exp_id, analyzers=[PfPRAnalyzer(working_dir=out_dir,
                                                                              sweep_variables=["Site_Name",
                                                                                               "Run_Number",
                                                                                               "x_Temporary_Larval_Habitat",
                                                                                               "ACT_Coverage",
                                                                                               "IRS_Coverage",
                                                                                               "ITN_Coverage"
                                                                                               ])])

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
