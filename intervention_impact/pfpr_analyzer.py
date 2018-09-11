import pandas as pd
import os
import pdb

from simtools.SetupParser import SetupParser
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Utilities.COMPSUtilities import exps_for_suite_id
from simtools.Utilities.COMPSUtilities import COMPS_login


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

    run_type = "exp"

    if run_type=="exp":
        exps = {"initial": {"initial_burnin_4":"8c066962-e6b5-e811-a2c0-c4346bcb7275"},
            "interactions/gates_examples": {"gates_reruns_bugfix": "5858cb46-f3b5-e811-a2c0-c4346bcb7275"
                                }}

        for subfolder, int_list in exps.items():
            for int_name, exp_id in int_list.items():
                print("exp id is " + exp_id)
                colname = "initial_prev" if "initial" in int_name else "final_prev"
                am = AnalyzeManager(exp_id,
                                    analyzers=PfPRAnalyzer(expname=int_name,
                                                           colname=colname,
                                                           outdir=os.path.join(out_dir, subfolder),
                                                           sweep_variables=["Site_Name",
                                                                            "Run_Number",
                                                                            "x_Temporary_Larval_Habitat",
                                                                            "ACT_Coverage",
                                                                            "IRS_Coverage",
                                                                            "ITN_Coverage"
                                                                            ]
                                                           ),
                                    force_analyze=True)
                print(am.experiments)
                am.analyze()
