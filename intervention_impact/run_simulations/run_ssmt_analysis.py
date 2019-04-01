"""
This script is executed as entrypoint in the docker SSMT worker.
Its role is to collect the experiment ids and analyzers (w/o paramenters) and use Assets files and run the SSMTAnalysis in comps.
"""
from simtools.Analysis.SSMTAnalysis import SSMTAnalysis
from simtools.SetupParser import SetupParser
from pfpr_analyzer_ssmt import PfPRAnalyzer
import pandas as pd
from simtools.AssetManager.FileList import FileList

# Add these 2 files to server under current working dir's Assets dir, these 2 files not used by server, just for record
asset_files = FileList(root='.', files_in_root=["run_ssmt_analysis.py"])
asset_files.add_file("site_details.csv")

if __name__ == "__main__":
    SetupParser.default_block = 'HPC'
    SetupParser.init()

    sites = pd.read_csv("site_details.csv")

    experiments = {"interactions": "1e4c19f9-e353-e911-a2c0-c4346bcb1554"}

    tags = tags = {'Demo': 'dtktools SSMTAnalysis', 'WorkItem type': 'Docker'}

    for dirname, exp_id in experiments.items():
        args = dict(dir_name=dirname, report_names=sites["name"].tolist(),
                    sweep_variables=["Run_Number", "x_Temporary_Larval_Habitat",
                                     "CM_Drug", "CM_Coverage",
                                     "IRS_Coverage",
                                     "ITN_Coverage", "ITN_Start",
                                     "ATSB_Coverage", "ATSB_Initial_Effect",
                                     "Larvicide_Coverage",
                                     "Ivermectin_Coverage", "Ivermectin_Duration", "Ivermectin_Monthly_Rounds",
                                     "PEV_Coverage", "PEV_Waning_Config_class", "PEV_Waning_Config_Decay_Time_Constant",
                                     "TBV_Coverage", "TBV_Waning_Config_Decay_Time_Constant"])
        # no need to set working_dir. using default working_dir as '.'
        analysis = SSMTAnalysis(experiment_ids=[exp_id], analyzers=[PfPRAnalyzer], analyzers_args=[args],
                                asset_files=asset_files, tags=tags)

        analysis.analyze()