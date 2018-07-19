import pandas as pd
import os

from simtools.SetupParser import SetupParser
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer


class PfPRAnalyzer(BaseAnalyzer):


    def __init__(self, expname, sweep_variables=["Run_Number"], outdir=".", colname=None):
        super(PfPRAnalyzer, self).__init__()
        self.sweep_variables = sweep_variables
        self.filenames = ["output/MalariaSummaryReport_AnnualAverage.json"]
        self.expname = expname
        self.colname = self.expname if not colname else colname
        self.outdir = outdir


    def select_simulation_data(self, data, simulation):
        channeldata = data[self.filenames[0]]["DataByTime"]["PfPR_2to10"]
        simdata = pd.DataFrame({self.colname: channeldata})
        simdata = simdata[-2:-1]

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata


    def finalize(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        d = pd.concat(selected).reset_index(drop=True)

        d.to_csv(os.path.join(self.outdir, "{name}.csv".format(name=self.expname)), index=False)


if __name__ == "__main__":
    SetupParser.init("HPC")
    out_dir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
                           'map_intervention_impact', 'lookup_tables')
    am = AnalyzeManager("4683d052-b58a-e811-a2c0-c4346bcb7275",
                        analyzers=PfPRAnalyzer(expname="initial_prev_better_spacing",
                                               colname="initial_prev",
                                               outdir = out_dir,
                                               sweep_variables=["Site_Name",
                                                                "Run_Number",
                                                                "x_Temporary_Larval_Habitat",
                                                                "ITN_Coverage",
                                                                "ITN_Distributions",
                                                                "outbreak_fraction"
                                                               
                                                                ]
                                               ))
    am.analyze()
