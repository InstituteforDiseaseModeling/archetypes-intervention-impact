import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import os
from plotting.colors import load_color_palette
mpl.rcParams['pdf.fonttype'] = 42

homedir = 'E:/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/corr_itn_irs'
datadir = homedir
plotdir = homedir


def pct5(x):
    return np.percentile(x, 5)


def pct95(x):
    return np.percentile(x, 95)


def clean_datafile(datadir, merged_fname) :

    fname1 = 'ITN_IRS_correlated.csv'
    fname2 = 'ITN_IRS_uncorrelated.csv'

    df = pd.read_csv(os.path.join(datadir, fname1))
    df2 = pd.read_csv(os.path.join(datadir, fname2))
    df = pd.concat([df, df2])

    colnames = [x for x in df.columns.values if 'Coverage' in x]
    df = df.fillna(0)

    mdf = df.groupby(colnames + ['x_Temporary_Larval_Habitat', 'Site_Name'])['final_prev'].agg([np.mean, pct5, pct95]).reset_index()
    mdf = mdf.rename(columns={'mean' : 'PfPR',
                              'pct5' : 'PfPR 5',
                              'pct95' : 'PfPR 95'})
    mdf = mdf.sort_values(by=['Site_Name', 'x_Temporary_Larval_Habitat'])
    mdf.to_csv(os.path.join(datadir, merged_fname), index=False)


def plot_traces(df):

    mchannel = 'PfPR'
    lchannel = 'PfPR 5'
    uchannel = 'PfPR 95'
    colnames = [x for x in df.columns.values if 'Coverage' in x]

    sns.set_style('whitegrid')
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.07, top=0.9)
    palette = load_color_palette()

    for s, (site, sdf) in enumerate(df.groupby('Site_Name')):
        baseline = sdf[(sdf['ITN_Coverage'] == 0) & (sdf['IRS_Coverage'] == 0) & (sdf['ITN_IRS_Coverage'] == 0)]
        ax = fig.add_subplot(3,3,s+1)
        for gn, (g, gdf) in enumerate(sdf.groupby(colnames)):
            ax.plot(baseline[mchannel], gdf[mchannel], color=palette[gn], label=g)
            ax.fill_between(baseline[mchannel], gdf[lchannel], gdf[uchannel],
                            linewidth=0, alpha=0.3, color=palette[gn])
        ax.set_title(site)
        ax.set_xlim(-0.05, 0.85)
        ax.set_ylim(-0.05, 0.85)
        if s > 4 :
            ax.set_xlabel('initial PfPR 2to10')
        if s%3 == 0 :
            ax.set_ylabel('final PfPR 2to10')
    ax.legend()
    fig.suptitle(colnames)
    fig.savefig(os.path.join(plotdir, 'itn_irs_correlated.png'))
    fig.savefig(os.path.join(plotdir, 'itn_irs_correlated.pdf'), format='PDF')
    plt.show()


def plot_cases_averted(data_fname, pfprdf):

    df = pd.read_csv(data_fname)
    colnames = [x for x in df.columns.values if 'Coverage' in x]
    channel = 'New_Clinical_Cases'
    df[channel] = df[channel]/df['Population']*1000/3
    df = df.groupby(['Site_Name', 'x_Temporary_Larval_Habitat'] + colnames)[channel].agg([np.mean, pct5, pct95]).reset_index()
    df = df.rename(columns={'mean' : channel,
                              'pct5' : '%s 5' % channel,
                              'pct95' : '%s 95' % channel})

    df = df.sort_values(by=['Site_Name', 'x_Temporary_Larval_Habitat'])

    sns.set_style('whitegrid')
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.07, top=0.9)
    palette = load_color_palette()

    int_order = ['itn', 'irs', 'uncorr', 'corr']

    for s, (site, sdf) in enumerate(df.groupby('Site_Name')):
        ax = fig.add_subplot(3,3,s+1)
        baseline = sdf[(sdf['ITN_Coverage'] == 0) & (sdf['IRS_Coverage'] == 0) & (sdf['ITN_IRS_Coverage'] == 0)]
        basepfpr = pfprdf[(pfprdf['ITN_Coverage'] == 0) & (pfprdf['IRS_Coverage'] == 0) &
                          (pfprdf['ITN_IRS_Coverage'] == 0) & (pfprdf['Site_Name'] == site)]

        with_nets = sdf[(sdf['IRS_Coverage'] == 0) & (sdf['ITN_Coverage'] == 0.6)]
        with_irs = sdf[(sdf['IRS_Coverage'] == 0.6) & (sdf['ITN_Coverage'] == 0)]
        with_uncorr = sdf[(sdf['IRS_Coverage'] == 0.6) & (sdf['ITN_Coverage'] == 0.6)]
        with_corr = sdf[(sdf['ITN_IRS_Coverage'] == 0.6)]

        for i, intdf in enumerate([with_nets, with_irs, with_uncorr, with_corr]) :
            diff = [x - y for x, y in zip(baseline[channel], intdf[channel])]
            ax.plot(basepfpr['PfPR'], diff, '-', color=palette[i], label=int_order[i])

        ax.set_title(site)
        if s%3 == 0 :
            ax.set_ylabel('%s averted' % channel)
        if s > 4 :
            ax.set_xlabel('initial PfPR2to10')
        ax.set_xlim(-0.05, 0.85)
        ax.set_ylim(-50, 1200)

    ax.legend()
    plt.savefig(os.path.join(plotdir, 'itn_irs_cases_averted_by_site.png'))
    plt.savefig(os.path.join(plotdir, 'itn_irs_cases_averted_by_site.pdf'), format='PDF')
    plt.show()


if __name__ == '__main__':

    # clean_datafile(datadir, 'itn_irs.csv')
    df = pd.read_csv(os.path.join(datadir, 'itn_irs.csv'))
    # plot_traces(df)
    plot_cases_averted(os.path.join(datadir, 'corr_itn_irs_cases.csv'), df)
