import os
import pandas as pd
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer

class PfPRAnalyzer(BaseAnalyzer):

    def __init__(self, dir_name, report_names=["AnnualAverage"], sweep_variables=None, working_dir=".",
                 last_year_only=True):
        super(PfPRAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.sitenames=report_names
        self.report_names = report_names  #Sharon added
        self.dir_name = dir_name
        self.last_year_only = last_year_only

    def select_simulation_data(self, data, simulation):
        colname = "initial_prev" if self.dir_name == "burnin" else "final_prev"

        simdata = []

        for site_name in self.sitenames:

            try:
                channeldata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTime"]["PfPR_2to10"]
                timedata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTime"]["Time Of Report"]
            except:
                raise FileNotFoundError("file not found for sim" + simulation.id)

            if self.last_year_only:
                tempdata = pd.DataFrame({colname: channeldata,
                                        "Site_Name": site_name})
                tempdata = tempdata[-2:-1]

            else:
                tempdata = pd.DataFrame({"day": timedata,
                                         colname: channeldata,
                                         "Site_Name": site_name})
                tempdata = tempdata[:-1]

            simdata.append(tempdata)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def filter(self, simulation):

        if "ClonedToRerun" in simulation.tags.keys():
            return False
        else:
            return True

    def finalize(self, all_data):
        data_sets_per_experiment = {}

        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.exp_name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)

            #sharon added following 2 lines to create dir with dir_name in comps server
            if not os.path.exists(os.path.join(self.working_dir, self.dir_name)):
                os.mkdir(os.path.join(self.working_dir, self.dir_name))

            d.to_csv(os.path.join(self.working_dir, self.dir_name, "{name}.csv".format(name=experiment_name)), index=False)


