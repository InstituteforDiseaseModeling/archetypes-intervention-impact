import os
import pdb

import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import exps_for_suite_id

from sim_output_processing.spatial_output_dataframe import construct_spatial_output_df

class BurdenAnalyzer(BaseAnalyzer):

    def __init__(self, dir_name, report_names=["AnnualAverage"], sweep_variables=None, working_dir="."):
        super(BurdenAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names] +
                                                  ['output/SpatialReport_%s.bin' % s for s in ["Population",
                                                                                        "Blood_Smear_Parasite_Prevalence",
                                                                                        'New_Infections',
                                                                                        'New_Clinical_Cases']]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.sitenames=report_names
        self.dir_name = dir_name
        self.reference = pd.read_csv('site_details.csv')
        self.spatial_channels = ["Population",
                                 "Blood_Smear_Parasite_Prevalence",
                                 'New_Infections',
                                 'New_Clinical_Cases']

    def select_simulation_data(self, data, simulation):
        simdata = []
        spatialdata = construct_spatial_output_df(data['output/SpatialReport_%s.bin' % self.spatial_channels[0]], self.spatial_channels[0])
        if len(self.spatial_channels) > 1:
            for ch in self.spatial_channels[1:]:
                spatialdata = pd.merge(left=spatialdata,
                                       right=construct_spatial_output_df(data['output/SpatialReport_%s.bin' % ch], ch),
                                       on=['time', 'node'])

        for site_name in self.sitenames:

            channeldata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTime"]["PfPR_2to10"]

            tempdata = pd.DataFrame({'PfPR2to10': channeldata,
                                    "Site_Name": site_name})
            tempdata = tempdata[-2:-1]

            nodeid = self.reference[self.reference['name'] == site_name]['nodeid'].values[0]
            for ch in self.spatial_channels :
                if ch in ['New_Infections', 'New_Clinical_Cases']:
                    tempdata[ch] = np.sum(spatialdata[spatialdata['node'] == nodeid][ch])
                else:
                    tempdata[ch] = np.mean(spatialdata[spatialdata['node'] == nodeid][ch])

            simdata.append(tempdata)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected).reset_index(drop=True)
        df.to_csv(os.path.join(self.working_dir, '%s.csv' % self.dir_name), index=False)


if __name__ == "__main__":
    SetupParser.init("HPC")
    out_dir = os.path.join('E:/', 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
                           'map_intervention_impact', 'lookup_tables')

    sites = pd.read_csv("site_details.csv")

    experiments = {
                   "corr_itn_irs" :"e7fd04a2-a8d8-e811-a2bd-c4346bcb1555"
                   }

    for dirname, exp_id in experiments.items():

        am = AnalyzeManager(exp_list=exp_id, analyzers=[BurdenAnalyzer(working_dir=out_dir,
                                                                     dir_name=dirname,
                                                                     report_names = sites["name"].tolist(),
                                                                      sweep_variables=["Run_Number",
                                                                                       "x_Temporary_Larval_Habitat",
                                                                                       "IRS_Coverage",
                                                                                       "ITN_Coverage",
                                                                                       "ITN_IRS_Coverage",
                                                                                       ])],
                            force_analyze=True)

        print(am.experiments)
        am.analyze()

