# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import colormaps
import matplotlib.cm as cm
import seaborn as sns
import collections as coll
import numpy as np
import scipy.stats as st
import os
import matplotlib.colors as cols
import functions as fxn

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def collapse_genecalls(recorded_call):
    """
    :param recorded_call: str of TCR V or J gene call(s), at the gene*allele level
    :return: list of unique gene calls (minus any orphon genes)
    """
    if ',' in recorded_call:
        call_list = [x for x in recorded_call.split(',') if '/OR' not in x]
        return list(set([x.split('*')[0] for x in call_list]))
    else:
        return [recorded_call.split('*')[0]]


def get_ci(vals, ci=0.95):
    """
    :param vals: list of numeric values
    :param ci: desired confidence interval (default = 95%)
    :return: confidence interval values
    """
    return st.t.interval(ci, len(vals)-1, loc=np.mean(vals), scale=st.sem(vals))


plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial'})
out_dir = fxn.get_plot_dir('gene-usage')

# Determine IMGT functionality (of prototypical genes at least)
dd, dd_path, fl_path, year, date, source, release = fxn.get_most_recent_genedb()

functionality = coll.defaultdict(str)
with open(fl_path, 'r') as in_file:
    for read, seq, qualnull in fxn.readfq(in_file):
        if '*01' in read and '|Homo sapiens|' in read:
            if '|TRA' in read or '|TRB' in read or '|TRG' in read or 'TRD' in read:
                bits = read.split('|')
                gene = bits[1].split('*')[0]
                prod = bits[3].replace('(', '').replace(')', '').replace('[', '').replace(']', '')
                functionality[gene] = prod
            

# Then read in the actual collapsed/autoDCR'd data, or use a pre-existing version if available
dat_check = [x for x in os.listdir(fxn.ref_dir) if '_mikelov-usage-data.tsv' in x]
if dat_check:
    dat_check.sort()
    usage_fl = dat_check[-1]
    print("Reading in usage data from previous run: " + usage_fl)
    usage_dat = pd.read_csv(os.path.join(fxn.ref_dir, usage_fl), sep='\t')

else:
    print("Reading in usage data from annotated repertoires...")
    repertoires = [x for x in os.listdir(fxn.mikelov_dir) if x.startswith('coll-') and
                   (x.endswith('.tsv.gz') or x.endswith('.tsv'))]
    repertoires.sort()

    gene_counts = {}
    usage_dat = []
    cnt_dat = []

    for fl in repertoires:
        donor = fl.split('-')[2].split('.')[0]
        print('\t', donor)
        tmp_gene_counts = {'TRAD': {'V': coll.Counter(), 'J': coll.Counter()},
                           'TRB': {'V': coll.Counter(), 'J': coll.Counter()}}

        dat = pd.read_csv(os.path.join(fxn.mikelov_dir, fl), sep='\t', index_col=False)

        for row in dat.index:
            row_dat = dat.loc[row]
            v_gene = collapse_genecalls(row_dat['v_call'])
            j_gene = collapse_genecalls(row_dat['j_call'])

            if len(v_gene) == 1 and len(j_gene) == 1:
                if 'TRB' in v_gene[0]:
                    locus = 'TRB'
                elif 'TRA' in v_gene[0] or 'TRD' in v_gene[0]:
                    locus = 'TRAD'

                else:
                    raise IOError('Unable to determine locus for ' + str(v_gene))

                tmp_gene_counts[locus]['V'][v_gene[0]] += int(row_dat['duplicate_count'])
                tmp_gene_counts[locus]['J'][j_gene[0]] += int(row_dat['duplicate_count'])

                cnt_dat.append([donor, locus, v_gene[0], j_gene[0],
                                len(row_dat['junction_aa']),
                                int(row_dat['duplicate_count'])])

            else:
                continue

        gene_counts[donor] = tmp_gene_counts

        for locus in ['TRAD', 'TRB']:
            for r in ['V', 'J']:
                for gene in tmp_gene_counts[locus][r]:
                    if 'TRDV' in gene:
                        gtype = 'D'
                    elif '/DV' in gene:
                        gtype = 'A/D'
                    elif 'TRA' in gene:
                        gtype = 'A'
                    elif 'TRB' in gene:
                        gtype = 'B'
                    else:
                        raise IOError('Unexpected gene type: ' + gene)
                    usage_dat.append([donor, locus, gene, r, gtype, tmp_gene_counts[locus][r][gene],
                                      tmp_gene_counts[locus][r][gene]/sum(tmp_gene_counts[locus][r].values())*100])

    cnt_dat = fxn.list_to_df(cnt_dat, ['donor', 'locus', 'v', 'j', 'cdr3len', '#'], False)
    usage_dat = fxn.list_to_df(usage_dat, ['donor', 'locus', 'gene', 'region', 'chain', '#', '%'], False)

    cnt_dat.to_csv(os.path.join(fxn.ref_dir, fxn.today() + '_mikelov-count-data.tsv.gz'), index=False, sep='\t')
    usage_dat.to_csv(os.path.join(fxn.ref_dir, fxn.today() + '_mikelov-usage-data.tsv'), index=False, sep='\t')

# TR[AD]V first
# Plot manual stripplot-style plot, with jittered dots per gene with their own CIs
trad_dat = usage_dat.loc[(usage_dat['region'] == 'V') & (usage_dat['chain'].isin(['A', 'A/D', 'D']))]

genes = list(set(trad_dat['gene']))
genes.sort()
x_dict = {'A': 1, 'A/D': 2, 'D': 3}
col_dict = {'A': 'tab:green', 'A/D': 'tab:purple', 'D': 'tab:cyan',
            'F': 'tab:cyan', 'ORF': 'tab:purple:', 'P': 'tab:red'}

plot_dat = []
for g in genes:
    g_dat = trad_dat.loc[trad_dat['gene'] == g]
    chain = g_dat.iloc[0]['chain']
    y = np.mean(g_dat['%'])
    x = x_dict[chain]
    n = len(g_dat)
    if n > 1:
        ci_lower, ci_upper = get_ci(g_dat['%'])
    else:
        ci_lower, ci_upper = '', ''
    f = functionality[g]
    plot_dat.append([g, chain, n, x, y, ci_lower, ci_upper, f])

plot_dat = fxn.list_to_df(plot_dat,
                          ['gene', 'chain', 'n', 'x', 'y', 'ci_lower', 'ci_upper', 'functionality'],
                          False)
plot_dat.to_csv(os.path.join(fxn.ref_dir, fxn.today() + '_mikelov-plot1.tsv'), index=False, sep='\t')

marker_dict = {'F': 'o', 'ORF': '^', 'P': 'v'}
edge_dict = {'F': cols.colorConverter.to_rgba('aqua', alpha=.5),
             'ORF': cols.colorConverter.to_rgba('purple', alpha=.5),
             'P': cols.colorConverter.to_rgba('red', alpha=.5)}

used_functionalities = list(set([functionality[x] for x in plot_dat['gene']]))
used_functionalities = [x for x in marker_dict if x in used_functionalities]  # cheeky sort

m_size = 8
x_jitter = 0.35

fig, ax = plt.subplots(figsize=(4, 5))
norm = Normalize(vmin=0, vmax=140)
colourmap = colormaps.get_cmap('crest')

for row in plot_dat.index:
    row_dat = plot_dat.loc[row]
    f = functionality[row_dat['gene']]
    colour = colourmap(norm(row_dat['n']))
    edgecol = edge_dict[f]
    x = np.random.uniform(row_dat['x']-x_jitter, row_dat['x']+x_jitter)
    if row_dat['ci_upper']:
        plt.plot([x, x], [row_dat['ci_lower'], row_dat['ci_upper']], c=colour, alpha=.6)
    plt.plot(x, row_dat['y'], c=colour,  marker=marker_dict[f], markersize=m_size, alpha=.8,
             markeredgecolor=edgecol)

for x in [1.5, 2.5]:
    plt.axvline(x=x, linewidth=1, color='gray', alpha=0.1)

sns.despine(top=True, right=True)
plt.xticks(ticks=[1, 2, 3], labels=['AV', 'AV/DV', 'DV'], minor=False)
plt.yscale('log')
plt.xlim(0.5, 3.5)
plt.ylim(0.00004, 20)
plt.ylabel('% gene usage')
plt.xlabel('TR chain(s)')
plt.savefig(os.path.join(out_dir, 'trad-usage.png'), dpi=300, bbox_inches='tight')
plt.close()

# And then plot the rest!
rest_dat = usage_dat.loc[((usage_dat['region'] == 'J') | (usage_dat['chain'].isin(['B']))) &
                         (usage_dat['chain'] != 'D')]

genes = list(set(rest_dat['gene']))
genes.sort()
x_dict = {'BV': 1, 'BJ': 2, 'AJ': 3}

plot_dat = []
for g in genes:
    g_dat = rest_dat.loc[rest_dat['gene'] == g]
    chain = g_dat.iloc[0]['chain'] + g_dat.iloc[0]['region']
    y = np.mean(g_dat['%'])
    x = x_dict[chain]
    n = len(g_dat)
    if n > 1:
        ci_lower, ci_upper = get_ci(g_dat['%'])
    else:
        ci_lower, ci_upper = '', ''
    f = functionality[g]
    plot_dat.append([g, chain, n, x, y, ci_lower, ci_upper, f])

plot_dat = fxn.list_to_df(plot_dat, ['gene', 'chain', 'n', 'x', 'y', 'ci_lower', 'ci_upper', 'functionality'], False)
plot_dat.to_csv(os.path.join(fxn.ref_dir, fxn.today() + '_mikelov-plot2.tsv'), index=False, sep='\t')


used_functionalities = list(set([functionality[x] for x in plot_dat['gene']]))
used_functionalities = [x for x in marker_dict if x in used_functionalities]  # cheeky sort

fig, ax = plt.subplots(figsize=(4, 5))
norm = Normalize(vmin=0, vmax=140)

for row in plot_dat.index:
    row_dat = plot_dat.loc[row]
    f = functionality[row_dat['gene']]
    colour = colourmap(norm(row_dat['n']))
    edgecol = edge_dict[f]
    x = np.random.uniform(row_dat['x']-x_jitter, row_dat['x']+x_jitter)
    if row_dat['ci_upper']:
        plt.plot([x, x], [row_dat['ci_lower'], row_dat['ci_upper']], c=colour, alpha=.6)
    plt.plot(x, row_dat['y'], c=colour,  marker=marker_dict[f], markersize=m_size, alpha=.8,
             markeredgecolor=edgecol)

for x in [1.5, 2.5]:
    plt.axvline(x=x, linewidth=1, color='gray', alpha=0.1)

sns.despine(top=True, right=True)
plt.xticks(ticks=[1, 2, 3], labels=['BV', 'BJ', 'AJ'], minor=False)
plt.yscale('log')
plt.xlim(0.5, 3.5)
plt.ylim(0.00004, 20)
plt.ylabel('% gene usage')
plt.xlabel('TR chain(s)')
plt.savefig(os.path.join(out_dir, 'rest-usage.png'), dpi=300, bbox_inches='tight')
plt.close()

# Legend plotting (for assembling plots downstream)
fig, ax = plt.subplots(figsize=(8, 3))
norm = Normalize(vmin=0, vmax=140)
plt.axis('off')

# Add legend bits
for m in marker_dict:
    if m in used_functionalities:
        plt.plot(-10, 0, c='black', marker=marker_dict[m], markersize=int(m_size * 1.3), label=m, lw=0,
                 markeredgecolor=edge_dict[m])

plt.xlim(0.5, 3.5)
plt.ylim(0.00004, 20)
plt.colorbar(cm.ScalarMappable(norm=norm, cmap=colourmap), ax=ax, label='# donors',
             location='top', ticks=[x for x in range(0, 160, 20)], aspect=50)
ax.legend(loc='lower center')
plt.savefig(os.path.join(out_dir, 'legend.png'), dpi=300, bbox_inches='tight', transparent=True)
plt.close()
