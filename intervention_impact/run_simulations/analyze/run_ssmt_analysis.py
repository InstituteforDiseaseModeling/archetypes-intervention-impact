"""
This script is executed as entrypoint in the docker SSMT worker.
Its role is to collect the experiment ids and analyzers (w/o paramenters) and use Assets files and run the SSMTAnalysis in comps.
"""
from simtools.Analysis.SSMTAnalysis import SSMTAnalysis
from simtools.SetupParser import SetupParser
from pfpr_analyzer_ssmt import PfPRAnalyzer
from incidence_analyzer_ssmt import IncAnalyzer
import pandas as pd
import os
import json
from simtools.AssetManager.FileList import FileList

version_name = "20200506_reextract_20191009_mega_era5_new_arch"
main_dir = os.path.join(os.path.expanduser("~"),
                            "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                            version_name, "input")
with open(os.path.join(main_dir, "input_params.json")) as f:
    instructions = json.loads(f.read())

# Add these 2 files to server under current working dir's Assets dir, these 2 files not used by server, just for record
asset_files = FileList(root='.', files_in_root=["run_ssmt_analysis.py"])
asset_files.add_file(os.path.join(main_dir, "site_details.csv"))

if __name__ == "__main__":
    SetupParser.default_block = 'HPC'
    SetupParser.init()

    sites = pd.read_csv(os.path.join(main_dir, "site_details.csv"))

    experiments = {"burnin": instructions["burnin_id"],
                   "intervention": instructions["intervention_id"]
                   }

    tags = {'Demo': 'dtktools SSMTAnalysis', 'WorkItem type': 'Docker'}

    for dirname, exp_ids in experiments.items():
        args = dict(dir_name=dirname, report_names=sites["id"].tolist(),
                    sweep_variables=["Run_Number", "x_Temporary_Larval_Habitat",
                                     "CM_Drug", "CM_Coverage",
                                     "IRS_Coverage",
                                     "ITN_Coverage", "ITN_Start",
                                     "ATSB_Coverage", "ATSB_Initial_Effect",
                                     "Larvicide_Coverage",
                                     "Ivermectin_Coverage", "Ivermectin_Duration", "Ivermectin_Monthly_Rounds",
                                     "PEV_Coverage", "PEV_Waning_Config_class", "PEV_Waning_Config_Decay_Time_Constant",
                                     "TBV_Coverage", "TBV_Waning_Config_Decay_Time_Constant",
                                     "SMC_Coverage", "SMC_Start", "SMC_Adherence", "SMC_Max_Age"],
                    last_year_only=False)

        exp_ids = exp_ids if isinstance(exp_ids, list) else [exp_ids]

        # no need to set working_dir. using default working_dir as '.'
        analysis = SSMTAnalysis(experiment_ids=exp_ids, analyzers=[IncAnalyzer], analyzers_args=[args],
                                asset_files=asset_files, tags=tags, analysis_name="Intervention Impact Analysis")

        analysis.analyze()