# -*- coding: utf-8 -*-

import os.path
import pyalluvial.alluvial as alluvial
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import functions as fxn

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


plt.rcParams.update({'font.sans-serif': 'Arial', 'font.size': 18})
dat_f = os.path.join(fxn.ref_dir, 'systematic-tcrseq-germline-reference-review.csv')
out_dir = fxn.get_plot_dir('review-alluvial')

df = pd.read_csv(dat_f, encoding='latin')
df['freq'] = 1

subset = ['Species', 'Loci', 'Sequence production', 'Analysis tool',
          'Germline source', 'Germline resource', 'Specific version?']

wide_df = df[subset + ['freq']].groupby(subset, as_index=False).sum()

plt.rcParams.update({'font.sans-serif': 'Arial', 'font.size': 8, 'font.weight': 'bold'})

plot = alluvial.plot(df=wide_df, xaxis_names=subset, y_name='freq',
                     alluvium='Specific version?',
                     ignore_continuity=True, figsize=(10, 5))
plt.savefig(os.path.join(out_dir, 'alluvial.pdf'), dpi=500, bbox_inches='tight')
plt.close()


total_setups = sum(wide_df['freq'])
wide_df['pc'] = wide_df['freq']/total_setups * 100

wide_df.to_csv(os.path.join(out_dir, 'processed_wideDF.tsv'), sep='\t', index=False)

# Create table of what percentage of each field each value makes
pcs = []
for field in subset:
    entries = list(set(wide_df[field]))
    entries.sort()
    for entry in entries:
        sub_df = wide_df.loc[wide_df[field] == entry]
        total = sum(sub_df['freq'])
        totalpc = sum(sub_df['pc'])
        pcs.append([field, entry, total, totalpc])

pcs = fxn.list_to_df(pcs, ['field', 'value', 'count', 'pc'], False)
pcs.to_csv(os.path.join(out_dir, 'percentages.tsv'), sep='\t', index=False)

# And a summary graph (for the graphical abstract)
plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial', 'font.weight': 'bold',
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Arial:italic', 'mathtext.rm': 'Arial',
                     'mathtext.bf': 'Arial:bold', 'mathtext.default': 'bf'})

summary = fxn.list_to_df([
           ['No GL info', pcs.loc[(pcs['field'] == 'Germline source') & (pcs['value'] == 'NP')].iloc[0]['pc']],
           ['Source', sum(pcs.loc[(pcs['field'] == 'Germline source') & (pcs['value'] != 'NP')]['pc'])],
           ['Resource', sum(pcs.loc[(pcs['field'] == 'Germline resource') & (pcs['value'] != 'NP')]['pc'])],
           ['Specific', pcs.loc[(pcs['field'] == 'Specific version?') & (pcs['value'] == 'Y')].iloc[0]['pc']]
           ], ['field', '% of experimental setups'], False)

sns.barplot(data=summary, x='field', y='% of experimental setups')
field_order = ['No GL info', 'Source', 'Resource', 'Specific']


g = sns.catplot(data=summary, x='field', y='% of experimental setups', kind='bar',
                hue='field', order=field_order, legend=False,
                palette={'No GL info': 'tab:red',
                         'Source': 'tab:orange',
                         'Resource': 'tab:olive',
                         'Specific': 'tab:green'})

sns.despine(top=True, right=True)
g.set_xticklabels(rotation=90)
g.set(ylim=(0, fxn.pc_rounding(max(summary['% of experimental setups']))))
g.set_axis_labels('', '% of published experimental setups')
out_path = os.path.join(out_dir, 'summary-barplot.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close('all')


