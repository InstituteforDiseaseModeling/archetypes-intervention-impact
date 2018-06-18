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

df_cols = ['Run_Number', 'x_Temporary_Larval_Habitat', 'ITN_Coverage', "Healthseek_Coverage",
           "IRS_Coverage"]

serialization_exp_ids = {'initial': "476808e8-6369-e811-a2c0-c4346bcb7275",
                         'final': "0bf64336-2e6b-e811-a2c0-c4346bcb7275"}

out_fname = "moine_climate_burnin.csv"
COMPS_login("https://comps.idmod.org")

full_prev_list = []

for run_type, serialization_exp_id in serialization_exp_ids.items():

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
        prev_df = prev_df[['PfPR_2to10']]

        for col in df_cols:
            try:
                prev_df[col] = df[col][x]
            except KeyError:
                continue

        prev_list.append(prev_df)

    prev_all = pd.concat(prev_list)
    prev_all['run_type'] = run_type

    full_prev_list.append(prev_all)

before_after = pd.concat(full_prev_list)
means = before_after.groupby(['run_type', 'x_Temporary_Larval_Habitat']).mean().drop('Run_Number', axis=1).reset_index()
means = pd.pivot_table(means, values='PfPR_2to10', index=['x_Temporary_Larval_Habitat'], columns=['run_type'])
means = means.query('initial>0 | final>0')
means = means[['initial', 'final']]
means.to_csv(os.path.join(out_dir, 'lookup_table_v2.csv'))

# prev_all.to_csv(os.path.join(out_dir, out_fname), index=False)


print(means)