import os
import pandas as pd
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer

class IncAnalyzer(BaseAnalyzer):

    def __init__(self, dir_name, report_names=["AnnualAverage"], sweep_variables=None, working_dir=".",
                 last_year_only=True):
        super(IncAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.sitenames=report_names
        self.report_names = report_names  #Sharon added
        self.dir_name = dir_name
        self.last_year_only = last_year_only

    def select_simulation_data(self, data, simulation):
        inc_colname = "initial_inc" if self.dir_name == "burnin" else "final_inc"
        severe_colname = "initial_severe_inc" if self.dir_name == "burnin" else "final_severe_inc"

        simdata = []

        for site_name in self.sitenames:

            try:
                timeinterval = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["Metadata"]["Reporting_Interval"]
                incdata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTimeAndAgeBins"]["Annual Clinical Incidence by Age Bin"]
                severedata =  data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTimeAndAgeBins"]["Annual Severe Incidence by Age Bin"]
                popdata =  data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTimeAndAgeBins"]["Average Population by Age Bin"]
                timedata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTime"]["Time Of Report"]
            except:
                raise FileNotFoundError("file not found for sim" + simulation.id)

            max_years = int(max(timedata) / timeinterval)
            incidence = [np.average(incdata[x], weights=popdata[x]) for x in range(0, max_years)]
            severe_incidence = [severedata[x][0] for x in range(0, max_years)]

            if self.last_year_only:
                tempdata = pd.DataFrame({inc_colname: incidence,
                                         severe_colname: severe_incidence,
                                        "Site_Name": site_name})
                tempdata = tempdata[-1]

            else:
                tempdata = pd.DataFrame({"day": timedata[0:max_years],
                                         inc_colname: incidence,
                                         severe_colname:severe_incidence,
                                         "Site_Name": site_name})

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

            d.to_csv(os.path.join(self.working_dir, self.dir_name, "{name}_incidence.csv".format(name=experiment_name)),
                     index=False)


