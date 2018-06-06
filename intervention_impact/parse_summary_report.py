import os
import pandas as pd
import sys
import json
import pdb
import re
import shutil

from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

out_dir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
                       'map_intervention_impact', 'prelim_itn_sweeps')

df_cols = ['Run_Number', 'x_Temporary_Larval_Habitat', 'add_ITN_age_season.coverage_all', "Healthseek_Coverage",
           "IRS_Coverage"]

serialization_exp_id = "61ba48ed-1469-e811-a2c0-c4346bcb7275"
out_fname = "moine_climate_burnin.csv"

COMPS_login("https://comps.idmod.org")
expt = ExperimentDataStore.get_most_recent_experiment(serialization_exp_id)

df = pd.DataFrame([x.tags for x in expt.simulations])
df['outpath'] = pd.Series([sim.get_path() for sim in expt.simulations])

prev_list = []

for x in df.index:
    print(x)
    report_dir = os.path.join(df['outpath'][x], 'output', 'MalariaSummaryReport_AnnualAverage.json')
    # report_dir = os.path.join(df['outpath'][x], 'output', 'VectorSpeciesReport.json')


    with open(report_dir) as f:
        report = json.loads(f.read())

    prev_df = pd.DataFrame(report['DataByTime'])
    prev_df = prev_df[-2:-1]
    for col in df_cols:
        try:
            prev_df[col] = df[col][x]
        except KeyError:
            continue

    prev_list.append(prev_df)

prev_all = pd.concat(prev_list)
prev_all.to_csv(os.path.join(out_dir, out_fname), index=False)

means = prev_all.groupby('x_Temporary_Larval_Habitat').mean()
print(means)