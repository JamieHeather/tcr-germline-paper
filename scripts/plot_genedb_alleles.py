# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import seaborn as sns
import functions as fxn
import argparse
import collections as coll
import numpy as np

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(description="plot_genedb_alleles v" + str(__version__) + '\n' +
                                                 ": plot the history of IMGT/GENE-DB alleles")

    # Input and output options
    parser.add_argument('-s', '--species', required=False, type=str, default='Homo sapiens',
                        help="Species list. Case sensitive, comma delimited, in quotes. "
                             "Default = \"Homo sapiens\".")

    parser.add_argument('-l', '--loci', required=False, type=str, default='TRA,TRB,TRG,TRD',
                        help="Loci list. Case sensitive, comma delimited, three-character identifiers in quotes. "
                             "Default = \"TRA,TRB,TRG,TRD\".")

    parser.add_argument('-r', '--regions', required=False, type=str, default='V,J',
                        help="Region list. Case sensitive, comma delimited, one-character identifiers in quotes. "
                             "Default = \"V,J\".")

    parser.add_argument('--version', action='version', version=__version__,
                        help="Print current version.")

    return parser.parse_args()


plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial', 'font.weight': 'bold',
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Arial:italic', 'mathtext.rm': 'Arial',
                     'mathtext.bf': 'Arial:bold', 'mathtext.default': 'bf'})

colours = ['gray', 'darkorange', 'tomato', 'green', 'olivedrab']
lstyles = [':', '--', '-.', '-']


def plot_locus_history(locus_df, used_dates, plot_path):
    """
    :param locus_df: pandas dataframe detailing the contents of a specific IMGT/GENE-DB locus
    :param used_dates: list of iso-format dates, of when GENE-DB was sampled
    :param plot_path: str detailing the path to save the plots in
    :return: nothing, makes the relevant plot in plot_path
    """
    y_height = 0.7 + (len(locus_df) / 13)
    y_labels = []
    fig, ax = plt.subplots(figsize=(5.5, y_height))

    y = 0
    for a_row in locus_df.index:
        row_dat = locus_df.loc[a_row]
        y_labels.append(row_dat['gene/allele'])
        # Plot lines first
        for ln in row_dat['lines']:
            plt.plot([ln[0], ln[1]], [y, y],
                     c=ln[2], ls=ln[3], alpha=ln[4])

        # Then markers
        for mr in row_dat['markers']:
            plt.plot(mr[0], y, c=mr[1], marker=mr[2], alpha=mr[3], markersize=mr[4])

        y += 1

    # And plot manual minor x ticks to show sampling dates
    for d in used_dates:
        plt.plot(d, -2, '^', c='purple', markersize=8)

    sns.despine(top=True, right=True)
    ax.set_yticks(range(len(y_labels)), y_labels, fontsize=6)
    ax.set_xticks([fxn.iso(str(x) + '-01-01') for x in range(2012, 2028, 2)],
                  labels=[str(x) for x in range(2012, 2028, 2)])  # , fontsize=6)
    plt.ylim(-2, len(y_labels) + 1)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close('all')


def earliest_novel_date(novel_row, dates_dict):
    """
    Determine the earliest date with which a filtered novel alleles was used
    :param novel_row: pandas Series, taken from the df of novel alleles detailing a specific allele
    :param dates_dict: dictionary of iso dates of when the different novel experiments were published
    :return: the earliest date that the given allele was reported (str, iso)
    """
    matches = novel_row[[x for x in novel_row.index if x.split('-')[0]
                         in dates_dict.keys() and x.endswith('-Donor-Count')]]
    match_nams = list(matches.index[[not np.isnan(x) for x in matches]])
    all_dates = [dates_dict[x.split('-')[0]] for x in match_nams]
    return min(all_dates)


if __name__ == "__main__":

    # Identify all the data directories, one per banked IMGT/GENE-DB release
    data_dirs = [x for x in os.listdir(fxn.release_dir)
                 if os.path.isdir(os.path.join(fxn.release_dir, x)) and x.startswith('20')]
    data_dirs.sort()

    search_str = 'fasta-nt-WithoutGaps-F+ORF+allP'

    input_args = vars(args())
    species_list = input_args['species'].split(',')
    species_list.sort()
    loci = input_args['loci'].split(',')
    regions = input_args['regions'].split(',')

    details = '_'.join(['-'.join([x.replace(' ', '') for x in species_list]),
                        '-'.join(loci),
                        '-'.join(regions)])
    out_dir = fxn.get_plot_dir(details)

    summary_dat = []
    summary_long = []

    # Then loop over all species to be examined...
    for species in species_list:

        # ... and then loop through all data directories, reading in the relevant data
        dat = []
        for dd in data_dirs:

            # Check file exists
            search = [x for x in os.listdir(os.path.join(fxn.release_dir, dd)) if search_str in x]
            if len(search) == 1:
                fl = search[0]
            else:
                print('Missing file in ' + dd)
                continue

            # Loop through
            year = dd[:4]
            date, source, release = dd.replace('/', '').split('_')

            approx_date = fxn.iso(date)
            with open(os.path.join(fxn.release_dir, dd, fl), 'r') as in_file:
                for header, seq, qualnull in fxn.readfq(in_file):
                    if species in header:
                        bits = header.rstrip().split('|')
                        if bits[1][:3] in loci and bits[1][3] in regions:
                            gene, allele = bits[1].split('*')
                            locus = gene[:3]
                            region = gene[3]
                            dat.append([date, approx_date, source, release, locus, gene, region, allele, seq.upper()])

        if not dat:
            raise IOError("No sequences detected - ensure correct input parameters (e.g. species name)")

        dat = fxn.list_to_df(dat, ['date', 'exact date', 'source', 'release',
                                   'locus', 'gene', 'region', 'allele', 'sequence'], False)

        print('Plotting ' + species + '!')
        dates = list(set(dat['exact date']))
        dates.sort()

        earliest = min(dat['exact date'])

        # Then loop through loci/genes/alleles and plot first occurrence/changes
        for locus in loci:
            for region in regions:

                plot_dat = []  # Running store of the plotting details for each allele in this locus/region

                y = 0
                lr_dat = dat.loc[(dat['locus'] == locus) & (dat['region'] == region)]
                genes = list(set(lr_dat['gene']))
                # genes.sort(key=len)
                genes.sort()

                lr_earliest = min(lr_dat['exact date'])

                for gene in genes:

                    g_dat = lr_dat.loc[dat['gene'] == gene]
                    g_dat.sort_values(by='date')
                    alleles = list(set(g_dat['allele']))
                    alleles.sort()

                    for allele in alleles:
                        # List of tuples giving: (x (timepoint), y (row), colour, marker, alpha, size)
                        marker_deets = []
                        line_deets = []  # Similarly: (x_start, x_end, y, colour, linestyle, alpha)

                        a_dat = g_dat.loc[g_dat['allele'] == allele]

                        # Loop through all timepoints after this allele's first appearance
                        # Record whether the sequence appears/changes/disappears
                        appeared = False
                        ref_date = a_dat.iloc[0]['exact date']
                        ref_seq = a_dat.iloc[0]['sequence']
                        discovered_dates = dates[dates.index(ref_date):]
                        latest_date = dates[0]
                        col_x = 0
                        ls_x = 0

                        for timepoint in discovered_dates:
                            match = a_dat.loc[a_dat['exact date'] == timepoint]

                            if len(match) > 1:
                                de_dup = match.drop_duplicates()
                                if len(de_dup) < len(match):
                                    print("Note: complete FASTA entry duplication identified in release " + release +
                                          " for allele " + gene + "*" + allele + " - continuing.")
                                    match = de_dup
                                else:
                                    print("Warning: different sequences with same ID observed in release " + release +
                                          " for allele " + gene + "*" + allele + " - ignoring this entry.")
                                    continue

                            if len(match) == 1 and not appeared:
                                if timepoint == earliest:
                                    marker_deets.append((timepoint, 'green', 'o', .5, 5))
                                else:
                                    marker_deets.append((timepoint, 'green', '^', .5, 5))
                                x_start = timepoint
                                appeared = True

                            if len(match) == 1 and appeared:
                                if match.iloc[0]['sequence'] != ref_seq:
                                    x_end = timepoint
                                    line_deets.append((x_start, x_end, colours[col_x], lstyles[ls_x], 0.7))
                                    marker_deets.append((timepoint, 'm', '*', .7, 6))

                                    col_x += 1
                                    ls_x += 1
                                    x_start = timepoint
                                    ref_seq = match.iloc[0]['sequence']

                                if timepoint > latest_date:
                                    latest_date = timepoint

                            if len(match) == 0 and appeared:
                                marker_deets.append((timepoint, 'r', 'v', .5, 5))
                                x_end = timepoint
                                line_deets.append((x_start, x_end, colours[col_x], lstyles[ls_x], 0.7))
                                appeared = False

                            if len(match) == 1 and dates.index(timepoint) + 1 == len(dates) and appeared:
                                marker_deets.append((timepoint, 'b', 's', .5, 5))

                                x_end = timepoint
                                line_deets.append((x_start, x_end, colours[col_x], lstyles[ls_x], 0.7))

                        gene_allele = a_dat.iloc[0]['gene'] + '*' + a_dat.iloc[0]['allele']

                        marker_deets = list(set(marker_deets))
                        line_deets = list(set(line_deets))

                        # Determine the status of this gene
                        if discovered_dates[0] == earliest:
                            appeared = 'present at start'
                        elif discovered_dates[0] == lr_earliest:
                            appeared = 'present at species start'
                        else:
                            appeared = 'added later'

                        if latest_date == dates[-1]:
                            presence = 'present'
                        else:
                            presence = 'removed'

                        if len(line_deets) == 1:
                            persistence = 'unchanged'
                        else:
                            persistence = 'changed'

                        plot_dat.append([gene_allele, y, marker_deets, line_deets, appeared, presence, persistence])

                        y += 1

                plot_dat = fxn.list_to_df(plot_dat, ['gene/allele', 'y', 'markers', 'lines',
                                                     'appeared', 'presence', 'persistence'], False)

                summary_dat.append([species, locus, region, locus+region, len(genes), len(plot_dat),
                                    len(plot_dat.loc[plot_dat['appeared'] == 'added later']),
                                    len(plot_dat.loc[plot_dat['presence'] == 'removed']),
                                    len(plot_dat.loc[plot_dat['persistence'] == 'changed'])
                                    ])

                for field in [('appeared', 'added later'), ('presence', 'removed'), ('persistence', 'changed')]:
                    summary_long.append([species, locus, region, locus+region, len(genes), len(plot_dat),
                                         field[1], len(plot_dat.loc[plot_dat[field[0]] == field[1]]),
                                         len(plot_dat.loc[plot_dat[field[0]] == field[1]])/len(plot_dat)*100])

                out_path = os.path.join(out_dir, species.replace(' ', '') + '_' + locus + region + '.png')
                plot_locus_history(plot_dat, dates, out_path)

                changed = plot_dat.loc[(plot_dat['appeared'] == 'added later') |
                                       (plot_dat['presence'] == 'removed') |
                                       (plot_dat['persistence'] == 'changed')]

                if not changed.empty:
                    out_path = os.path.join(out_dir, species.replace(' ', '') + '_' + locus + region + '-ALTERED.png')
                    plot_locus_history(changed, dates, out_path)

    summary_dat = fxn.list_to_df(summary_dat, ['species', 'chain', 'region', 'locus',
                                               '# genes', '# alleles',
                                               '# present from start',
                                               '# present at current',
                                               '# unchanged'], False)

    summary_long = fxn.list_to_df(summary_long, ['species', 'chain', 'region', 'locus',
                                                 '# genes', '# alleles', 'field', 'count', '%'], False)
    field_order = ['added later', 'changed', 'removed']

    g = sns.catplot(data=summary_long, x='chain', y='%', col='field', kind='bar', row='species', hue='region',
                    order=loci, col_order=field_order)
    g.set_xticklabels(rotation=90)
    g.fig.subplots_adjust(hspace=0.1, wspace=0.1)
    g.set_titles("{col_name}")
    g.set(ylim=(0, fxn.pc_rounding(max(summary_long['%']))))
    g.set_axis_labels("", "% of total alleles")

    out_path = os.path.join(out_dir, 'summary.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close('all')

    for sp in species_list:
        subset = summary_long.loc[summary_long['species'] == sp]
        g = sns.catplot(data=subset, x='chain', y='%', col='field', kind='bar', hue='region', order=loci,
                        col_order=field_order, height=3, aspect=.55)
        g.set_xticklabels(rotation=90)
        g.fig.subplots_adjust(hspace=0.1, wspace=0.1)
        g.set_titles("{col_name}")
        g.set(ylim=(0, fxn.pc_rounding(max(subset['%']))))
        g.set_axis_labels("", "% of total alleles")
        out_path = os.path.join(out_dir, sp.replace(' ', '') + '-summary.png')
        plt.savefig(out_path, dpi=300, bbox_inches='tight')

        g = sns.catplot(data=subset, x='chain', y='%', col='field', kind='bar', hue='region', order=loci,
                        col_order=field_order, height=3, aspect=.65, legend=False)
        g.set_xticklabels(rotation=90)
        g.fig.subplots_adjust(hspace=0.1, wspace=0.1)
        g.set_titles("{col_name}")
        g.set(ylim=(0, fxn.pc_rounding(max(subset['%']))))
        g.set_axis_labels("", "% of total alleles")
        out_path = os.path.join(out_dir, sp.replace(' ', '') + '-nolegend-summary.png')
        plt.savefig(out_path, dpi=300, bbox_inches='tight')

        plt.close('all')

        # Single summary plot of just 'indeterminate alleles' per locus
        # (i.e. those which undergo some change in the tracked releases)
        summed_subset = []
        for locus in loci:
            for region in ['V', 'J']:
                specific_subset = subset.loc[(subset['chain'] == locus) & (subset['region'] == region)]
                summed_subset.append([locus, region,
                                      sum(specific_subset['count']) / specific_subset.iloc[0]['# alleles'] * 100])

        summed_subset = fxn.list_to_df(summed_subset, ['chain', 'region', '%'], False)

        g = sns.catplot(data=summed_subset, x='chain', y='%', kind='bar', hue='region', order=loci, legend=False)
        g.set_xticklabels(rotation=90)
        g.set_axis_labels("", "% alleles indeterminate")
        out_path = os.path.join(out_dir, sp.replace(' ', '') + '-total-indeterminate-summary.png')
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close('all')

# Plot a longitudinal summary of the number of alleles and genes, supplemented with novel alleles from other sources

# Need to establish the dates at which different novel alleles get discovered
# Only counting those found in 2 studies (but count from first date of discovery)
novel = fxn.get_novel_alleles(1, 1, '')

novel_dates = coll.defaultdict(str)
novel_dates_path = os.path.join(fxn.ref_dir, 'novel-publication-dates.tsv')

with open(novel_dates_path, 'r') as in_file:
    for line in in_file:
        bits = line.rstrip().split('\t')
        if len(bits) == 2:
            novel_dates[bits[0].split('-')[0]] = fxn.iso(bits[1])


novel_appearances = {}
novel_dates_sorted = {}
for locus in loci:
    novel_appearances[locus] = coll.Counter()
    locus_dat = novel.loc[novel['Gene'].str.startswith(locus)]
    for row in locus_dat.index:
        row_dat = locus_dat.loc[row]
        appearance = earliest_novel_date(row_dat, novel_dates)
        novel_appearances[locus][appearance] += 1
    novel_dates_sorted[locus] = list(novel_appearances[locus].keys())
    novel_dates_sorted[locus].sort()

sum_dat = []
novel_counts = coll.Counter()
for date in dates:
    date_dat = dat.loc[dat['exact date'] == date]
    release = date_dat.iloc[0]['release']
    for locus in loci:

        if locus not in novel_counts:
            novel_counts[locus] = 0

        lr_dat = date_dat.loc[(date_dat['locus'] == locus)]
        gene_count = len(list(set(lr_dat['gene'])))
        allele_count = len(lr_dat)

        if novel_dates_sorted[locus]:
            if date >= novel_dates_sorted[locus][0]:
                novel_counts[locus] += novel_appearances[locus][novel_dates_sorted[locus][0]]
                novel_dates_sorted[locus].pop(novel_dates_sorted[locus].index(novel_dates_sorted[locus][0]))

        sum_dat.append([date, release, locus, gene_count, 'gene'])
        sum_dat.append([date, release, locus, allele_count, 'allele'])

        if novel_counts[locus] > 0:
            sum_dat.append([date, release, locus, allele_count + novel_counts[locus], 'allele+novel'])

sum_dat = fxn.list_to_df(sum_dat, ['date', 'release', 'locus', 'count', 'type'], False)


fig, ax = plt.subplots(figsize=(5.5, 3))
for typ in ['gene', 'allele', 'allele+novel']:
    for chain in loci:
        subset = sum_dat.loc[(sum_dat['type'] == typ) & (sum_dat['locus'] == chain)]
        for idx in range(len(subset)):
            if idx + 1 < len(subset):
                x1 = subset.iloc[idx]['date']
                x2 = subset.iloc[idx+1]['date']
                y1 = subset.iloc[idx]['count']
                y2 = subset.iloc[idx+1]['count']
                if y2 > y1:
                    col = 'blue'
                elif y2 == y1:
                    col = 'black'
                else:
                    col = 'red'
                ax.plot([x1, x2], [y1, y2], c=col, alpha=.5, ls='--')


sns.scatterplot(data=sum_dat, x='date', y='count', style='locus', hue='type',
                s=100, linewidth=0, alpha=0.5, ax=ax, clip_on=False)
ax.set_ylabel('Count', weight='bold')
ax.set_xlabel('Year', weight='bold')
ax.set_xticks([fxn.iso(str(x)+'-01-01') for x in range(2012, 2028, 2)],
              labels=[str(x) for x in range(2012, 2028, 2)])
plt.ylim(50, 300)
sns.despine(top=True, right=True)
plt.legend(bbox_to_anchor=(1.07, 1), loc=2, borderaxespad=0., markerscale=2)
out_path = os.path.join(out_dir, '_'.join([species.replace(' ', '-'), ''.join([x[2] for x in loci])])) + '.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close('all')

ms = 12
fig, ax = plt.subplots(figsize=(1, 1))
ax.plot([-1], [-1], 'go', alpha=.5, markersize=ms, label='Present in earliest')
ax.plot([-1], [-1], 'g^', alpha=.5, markersize=ms, label='Appears later')
ax.plot([-1], [-1], 'm*', alpha=.7, markersize=ms, label='Sequence changes')
ax.plot([-1], [-1], 'rv', alpha=.5, markersize=ms, label='Allele disappears')
ax.plot([-1], [-1], 'bs', alpha=.5, markersize=ms, label='Still present')
ax.plot([-1], [-1], '^', c='purple', markersize=ms, label='GENE-DB sampled')
plt.xlim(10, 20)
plt.ylim(10, 20)
sns.despine(top=True, bottom=True, left=True, right=True)
ax.legend(loc=2, borderaxespad=0, ncols=2,
          columnspacing=0.1, handletextpad=0)
plt.axis('off')
out_path = os.path.join(out_dir, 'timeline-legend.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight', transparent=True)
plt.close()

fig, ax = plt.subplots(figsize=(1, 1))
ax.plot([-1], [-1], c='tab:blue', marker='s', alpha=.9, markersize=ms, label='gene', lw=0)
ax.plot([-1], [-1], c='tab:orange', marker='s', alpha=.9, markersize=ms, label='allele', lw=0)
ax.plot([-1], [-1], c='tab:green', marker='s', alpha=.9, markersize=ms, label='allele+novel', lw=0)
ax.plot([-1], [-1], c='black', marker='o', alpha=.9, markersize=ms, label='TRA', lw=0)
ax.plot([-1], [-1], c='black', marker='X', alpha=.9, markersize=ms, label='TRB', lw=0)
plt.xlim(10, 20)
plt.ylim(10, 20)
sns.despine(top=True, bottom=True, left=True, right=True)
ax.legend(loc=2, borderaxespad=0, ncols=2,
          columnspacing=0.3, handletextpad=0.01)
plt.axis('off')
out_path = os.path.join(out_dir, 'lineplot-legend.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight', transparent=True)
plt.close()


fig, ax = plt.subplots(figsize=(1, 1))
ax.plot([-1], [-1], c='tab:blue', marker='s', alpha=.9, markersize=ms, label='V', lw=0)
ax.plot([-1], [-1], c='tab:orange', marker='s', alpha=.9, markersize=ms, label='J', lw=0)
plt.xlim(10, 20)
plt.ylim(10, 20)
sns.despine(top=True, bottom=True, left=True, right=True)
ax.legend(loc=2, handletextpad=0.001, handlelength=1)
plt.axis('off')
out_path = os.path.join(out_dir, 'summary-legend.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight', transparent=True)
plt.close()
