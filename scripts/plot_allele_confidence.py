# -*- coding: utf-8 -*-

import os
import collections as coll
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import functions as fxn
import subprocess
from matplotlib.patches import Rectangle


__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial', 'font.weight': 'bold',
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Arial:italic', 'mathtext.rm': 'Arial',
                     'mathtext.bf': 'Arial:bold', 'mathtext.default': 'bf'})

colours = ['gray', 'darkorange', 'tomato', 'green', 'olivedrab']
lstyles = [':', '--', '-.', '-']


if __name__ == "__main__":

    data_dirs = [x for x in os.listdir(fxn.release_dir)
                 if os.path.isdir(os.path.join(fxn.release_dir, x)) and x.startswith('20')]
    data_dirs.sort()

    search_str = 'fasta-nt-WithoutGaps-F+ORF+allP'

    loci = ['TRA', 'TRB', 'TRG', 'TRD']
    locus_conv = {'TRA': 'TRAD', 'TRB': 'TRB', 'TRG': 'TRG', 'TRD': 'TRAD'}
    regions = ['V']

    out_dir = fxn.get_plot_dir('allele-confidence')

    species = 'Homo sapiens'
    dd, dd_path, fl_path, year, date, source, release = fxn.get_most_recent_genedb()

    approx_date = fxn.iso(date)
    tmp_hold = {}

    # Read in the GENE-DB data, and filter out the human genes for the requested loci/regions
    out_str = ''
    with open(fl_path, 'r') as in_file:
        for header, seq, qualnull in fxn.readfq(in_file):
            if species in header:
                bits = header.rstrip().split('|')

                if bits[1][:3] in loci and bits[1][3] in regions:
                    gene, allele = bits[1].split('*')
                    locus = gene[:3]
                    region = gene[3]

                    # Just hold on to the out row for now, as we also need to get the gapped sequence too
                    tmp_hold[bits[1]] = ['IMGT', date, approx_date, source, release, locus, locus_conv[locus],
                                         gene, region, allele, bits[1], seq.upper()]
                    out_str += fxn.fastafy(bits[1], seq.upper())

    # Write that out
    imgt_file = 'IMGT-human-ungapped.fasta'
    with open(imgt_file, 'w') as outf:
        outf.write(out_str)

    # Then make a copy and add novel alleles
    imgt_plusnovel_file = 'IMGTplusnovel-human-ungapped.fasta'
    shutil.copyfile(imgt_file, imgt_plusnovel_file)
    fxn.get_novel_alleles(2, 2, imgt_plusnovel_file)

    # Then generate gapped versions of each of those
    # Generate the individual locus-specific references, then combine into one and use as the reference
    for locus in loci:
        cmd = ('extract_refs -r ' +
               os.path.join(dd_path, 'IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP') +
               ' -L ' + locus + ' -F "Homo sapiens"')
        subprocess.call(cmd, shell=True)

    subprocess.call("ls Homo_sapiens_TR*V_gapped.fasta ", shell=True)
    subprocess.call("cat Homo_sapiens_TR*V_gapped.fasta > combined_TRxV_gapped.fasta", shell=True)

    files = [imgt_file, imgt_plusnovel_file]
    for fl in files:
        cmd = "gap_sequences " + fl + " combined_TRxV_gapped.fasta " + fl.replace('ungapped', 'gapped')
        subprocess.call(cmd, shell=True)

    subprocess.call("rm Homo*fasta", shell=True)

    # Having generated the relevant files, read the data in
    data = {'IMGT': {'ungapped': coll.defaultdict(str), 'gapped': coll.defaultdict(str)},
            'IMGTplusnovel': {'ungapped': coll.defaultdict(str), 'gapped': coll.defaultdict(str)}}

    for fl in files + [x.replace('ungapped', 'gapped') for x in files]:
        nam = fl.split('-')[0]
        typ = fl.split('-')[-1].split('.')[0]
        with open(fl, 'r') as in_file:
            for header, read, null in fxn.readfq(in_file):
                if header[3] == 'V':
                    data[nam][typ][header] = read

    # Trim all data back to CYS-104 (leaving that residue's codon in place)
    filt = {'IMGT': {'ungapped': coll.defaultdict(str), 'gapped': coll.defaultdict(str)},
            'IMGTplusnovel': {'ungapped': coll.defaultdict(str), 'gapped': coll.defaultdict(str)}}
    vlen = 312
    for typ in data:
        for ga in data[typ]['gapped']:
            to_trim = data[typ]['gapped'][ga][fxn.aa2nt(105):]
            if to_trim:
                if '.' in to_trim:
                    raise IOError("Gap detected in to_trim sequence! " + to_trim)

                if data[typ]['ungapped'][ga].endswith(to_trim):
                    trim = -len(to_trim)
                    seq = data[typ]['gapped'][ga][:trim]

                    # Filter out any short Vs
                    if len(seq) == vlen:
                        for gap in ['ungapped', 'gapped']:
                            filt[typ][gap][ga] = data[typ][gap][ga][:trim]
                    else:
                        print(ga + ' too short!')

    # Then go through each locus/allele, and ask for each 3'-fixed substr, how many gene/alleles/proteins does it match
    region_order = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3']
    dat = []
    long = []
    translation_counts = coll.Counter()
    for typ in ['IMGT', 'IMGTplusnovel']:
        for locus in loci:
            all_alleles = [x for x in list(filt[typ][gap].keys()) if locus in x or '/' + locus[-1]+'V' in x]
            ungapped = dict(zip(all_alleles, [filt[typ]['ungapped'][x] for x in all_alleles]))
            gapped = dict(zip(all_alleles, [filt[typ]['gapped'][x] for x in all_alleles]))
            gap_seqs = [gapped[x] for x in all_alleles]
            ungap_seqs = [ungapped[x] for x in all_alleles]
            ungap_seqs_nr = list(set(ungap_seqs))  # non redundant
            all_genes = list(set([x.split('*')[0] for x in gapped.keys()]))
            # translated = dict(zip(all_alleles, [translate(filt[typ]['ungapped'][x]) for x in all_alleles]))
            all_translations = [fxn.translate(filt[typ]['ungapped'][x]) for x in list(set(all_alleles))]
            translation_counts[locus + '-total'] = len(all_translations)
            translation_counts[locus + '-uniq'] = len(list(set(all_translations)))

            for al in all_alleles:
                gene, allele = al.split('*')
                # translations = [translate(filt[typ]['ungapped'][x]) for x in list(set(all_alleles))]

                for start in [x for x in range(0, vlen)][::-1]:

                    substr = gapped[al][start:].replace('.', '')
                    shared = [typ, locus, al, gene, allele, start, fxn.nt2region(start)]

                    # Allele ID level matching
                    matched_alleles = [x for x in ungap_seqs if substr in x]
                    matched_alleles_num = len(matched_alleles)
                    matched_alleles_pc = matched_alleles_num / len(ungap_seqs) * 100
                    allele_conf_pc = 1 / matched_alleles_num * 100
                    long.append(shared + ['allele (nt)', matched_alleles_num, matched_alleles_pc, allele_conf_pc])

                    # Gene level matching
                    matched_genes = list(set([x.split('*')[0] for x in ungapped if ungapped[x] in matched_alleles]))
                    matched_genes_num = len(matched_genes)
                    matched_genes_pc = matched_genes_num / len(all_genes) * 100
                    gene_conf_pc = 1 / matched_genes_num * 100
                    long.append(shared + ['gene', matched_genes_num, matched_genes_pc, gene_conf_pc])

                    # Translation level matching
                    translations = [fxn.translate(x) for x in list(set(matched_alleles))]
                    matched_aa_num = len(translations)
                    matched_aa_pc = matched_aa_num / len(all_translations) * 100
                    aa_conf_pc = 1 / matched_aa_num * 100
                    long.append(shared + ['allele (aa)', matched_aa_num, matched_aa_pc, aa_conf_pc])

                    dat.append(shared +
                               [matched_alleles_num, matched_alleles_pc, allele_conf_pc,
                                matched_genes_num, matched_genes_pc, gene_conf_pc,
                                matched_aa_num, matched_aa_pc, aa_conf_pc])

    dat = fxn.list_to_df(dat, ['data', 'locus', 'gene*allele', 'gene', 'allele', 'position', 'region',
                               '# allele matches', '% allele matches', 'allele confidence %',
                               '# gene matches', '% gene matches', 'gene confidence %',
                               '# translation matches', '% translation matches', 'translation confidence %'], False)

    long = fxn.list_to_df(long,
                          ['data', 'locus', 'gene*allele', 'gene', 'allele', 'position', 'region',
                           'sequence type', 'matched #', 'matched %', 'confidence %'], False)

    dat.to_csv(os.path.join(out_dir, 'wide-df.tsv'), sep='\t', index=False)
    dat.to_csv(os.path.join(out_dir, 'long-df.tsv'), sep='\t', index=False)

    for typ in ['IMGT', 'IMGTplusnovel']:
        for locus in loci:

            seq_order = ['allele (nt)', 'gene', 'allele (aa)']

            tdat = long.loc[(long['data'] == typ) & (long['locus'] == locus)]

            fig, ax = plt.subplots(figsize=(5.5, 5))
            sns.lineplot(data=tdat, x='position', y='confidence %', ax=ax,
                         hue='sequence type',
                         hue_order=seq_order, linewidth=3, linestyle='--')

            lw = 6
            h = 105
            v_bits = [(fxn.aa2nt(1), fxn.aa2nt(26), 'FR1', 'tab:cyan'),
                      (fxn.aa2nt(27), fxn.aa2nt(38), 'CDR1', 'tab:purple'),
                      (fxn.aa2nt(39), fxn.aa2nt(55), 'FR2', 'tab:cyan'),
                      (fxn.aa2nt(56), fxn.aa2nt(65), 'CDR2', 'tab:purple'),
                      (fxn.aa2nt(66), fxn.aa2nt(103), 'FR3', 'tab:cyan'),
                      (fxn.aa2nt(103), fxn.aa2nt(105), '3', 'tab:purple')]

            for vb in v_bits:
                plt.gca().add_patch(Rectangle(
                    (vb[0], h), (vb[1] - vb[0]) + 3, lw, linewidth=1, edgecolor='black', facecolor=vb[3]))
                ax.text(fxn.halfway(vb[0], vb[1]) + 1.5, h + lw / 2, vb[2], size=10, c='white',
                        horizontalalignment='center', verticalalignment='center')

            ax.fill(0, 0, 'white')  # No idea why, but the rectangles don't plot without this
            ax.set_yticks([0, 20, 40, 60, 80, 100])
            ax.spines['left'].set_bounds(0, 100)
            sns.despine(top=True, right=True)
            plt.xlabel("sequence 5' start site")
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(os.path.join(out_dir,
                                     '-'.join(['combALL', 'confidence', locus, typ, 'lineplot', '.png'])),
                        dpi=300, bbox_inches='tight')
            plt.close()


for locus in loci:
    print(locus, "has", str(translation_counts[locus+'-total']), "total translated sequences, of which",
          str(translation_counts[locus+'-uniq']), "are unique (" +
          str(round(translation_counts[locus+'-uniq']/translation_counts[locus+'-total'] * 100, 1)) + "%)")
