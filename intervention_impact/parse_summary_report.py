import os
import pandas as pd
import sys
import json
import pdb
import re
import shutil

from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.Utilities.COMPSUtilities import COMPS_login

serialization_exp_id = "31e731e5-be64-e811-a2c0-c4346bcb7275"

COMPS_login("https://comps.idmod.org")
expt = ExperimentDataStore.get_most_recent_experiment(serialization_exp_id)

df = pd.DataFrame([x.tags for x in expt.simulations])
df['outpath'] = pd.Series([sim.get_path() for sim in expt.simulations])

prev_list = []

for x in df.index:
    print(x)
    report_dir = os.path.join(df['outpath'][x], 'output', 'MalariaSummaryReport_AnnualAverage.json')

    with open(report_dir) as f:
        report = json.loads(f.read())

    prev_df = pd.DataFrame(report['DataByTime'])
    prev_df['x_temp'] = df['x_Temporary_Larval_Habitat'][x]

    prev_list.append(prev_df)

prev_all = pd.concat(prev_list)
last_year = prev_all.loc[49]
means = last_year.groupby('x_temp').mean()