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

df_cols = ['Run_Number', 'x_Temporary_Larval_Habitat', 'funestus.Anthropophily', 'funestus.Indoor_Feeding_Fraction',
           'ITN_Coverage', "ACT_Coverage", "IRS_Coverage"]

parse_type="simple"

serialization_exp_ids = {'initial': "0911f27b-7877-e811-a2c0-c4346bcb7275",
                         'final': "c42fbd67-ba77-e811-a2c0-c4346bcb7275"
                         }

out_fname = 'lookup_table_bajonapo_multi_int.csv'
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

        if parse_type=="all":
            datasets  = {key: pd.DataFrame(value) for (key, value) in report['DataByTimeAndAgeBins'].items()}

            for name, dataset in datasets.items():
                new_name = name.replace(" ", "_")
                dataset.columns=report['Metadata']['Age Bins']
                dataset['year'] = range(len(dataset))
                datasets[name] = pd.melt(dataset, id_vars="year", value_name=new_name, var_name="age")

            prev_df = pd.merge(datasets['PfPR by Age Bin'], datasets['Average Population by Age Bin'])

        else:
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
    # prev_all['run_type'] = run_type
    prev_all.rename(columns={'PfPR_2to10':run_type}, inplace=True)

    full_prev_list.append(prev_all)

if parse_type=="all":
    max_burnin_year = max(full_prev_list[0]['year'])
    full_prev_list[1]['year'] = full_prev_list[1]['year'] + max_burnin_year
    full_prev_list[0] = full_prev_list[0].query('year<@max_burnin_year')
    before_after = pd.concat(full_prev_list)

else:
    before_after = pd.merge(full_prev_list[0], full_prev_list[1])

# before_after.to_csv(os.path.join(out_dir, 'full_moine_data.csv'), index=False)

means = before_after.groupby(['x_Temporary_Larval_Habitat',
                              'IRS_Coverage', 'ITN_Coverage', 'ACT_Coverage']).mean().drop('Run_Number', axis=1).reset_index()
# means = pd.pivot_table(means, values='PfPR_2to10', index=['x_Temporary_Larval_Habitat', 'funestus.Anthropophily'], columns=['run_type'])
means = means.query('initial>0 | final>0')
# means = means[['initial', 'final']]
means.to_csv(os.path.join(out_dir, out_fname), index=False)

print(means)