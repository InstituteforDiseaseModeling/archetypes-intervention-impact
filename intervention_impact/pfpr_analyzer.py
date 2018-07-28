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

    exps = {"initial": {"initial_prev": "3dc6771b-0291-e811-a2c0-c4346bcb7275"},
            "itn_only": {"itn_baseline": "49e54528-1891-e811-a2c0-c4346bcb7275",
                         "itn_corr_nets": "df13a542-1e91-e811-a2c0-c4346bcb7275"},
            "irs_only": {"irs_baseline": "aebafa27-2291-e811-a2c0-c4346bcb7275"},
            "act_only": {"act_baseline": "a69da6dd-2491-e811-a2c0-c4346bcb7275",
                         "act_hs_rate": "b9ea054c-2791-e811-a2c0-c4346bcb7275"}
            }

    for subfolder, int_list in exps.items():
        colname = "initial_prev" if subfolder=="initial" else "final_prev"
        for int_name, exp_id in int_list.items():

            am = AnalyzeManager(exp_id,
                                analyzers=PfPRAnalyzer(expname=int_name,
                                                       colname=colname,
                                                       outdir = os.path.join(out_dir, subfolder),
                                                       sweep_variables=["Site_Name",
                                                                        "Run_Number",
                                                                        "x_Temporary_Larval_Habitat",
                                                                        "ACT_Coverage",
                                                                        "IRS_Coverage",
                                                                        "ITN_Coverage",
                                                                        "Hates_Nets",
                                                                        "ACT_HS_Rate"
                                                                        ]
                                                       ),
                                force_analyze=True)
            am.analyze()
