"""
Do a simple collapse, purely based on exact matching V/J/CDR3
"""


import gzip
import pandas
import argparse
import collections as coll
import pandas as pd
import functions as fxn

def args():
    """
    :return: argparse details of command line arguments passed to the script
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description='')

    # Add arguments
    parser.add_argument('-in', '--in_file', type=str, required=True,
                        help='autoDCR output TSV file')

    parser.add_argument('-t', '--threshold', type=int, required=False, default=1,
                        help='Minimum count frequency to keep TCRs')

    return parser.parse_args()


def opener(file_path, open_mode):
    """
    :rtype: object
    :param file_path: path to file to be opened
    :param open_mode: mode by which to open the file (e.g. w/r/a)
    :return: the appropriate file opening command (open or gzip.open)
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, open_mode + 't')
    else:
        return open(file_path, open_mode)


out_headers = ['sequence_id', 'v_call', 'd_call', 'j_call', 'junction_aa', 'duplicate_count', 'sequence',
               'junction', 'rev_comp', 'productive', 'sequence_aa',
               'inferred_full_nt', 'inferred_full_aa', 'non_productive_junction_aa',
               'vj_in_frame', 'stop_codon', 'conserved_c', 'conserved_f',
               'inter_tag_seq', 'inter_tag_qual', 'umi_seq', 'umi_qual',
               'sequence_alignment', 'germline_alignment', 'v_cigar', 'd_cigar', 'j_cigar']


if __name__ == '__main__':

    # Determine the requested input/output parameters
    input_args = vars(args())

    dat = pd.read_csv(input_args['in_file'], sep='\t')
    # Get v/j/cdr3 info for productive rearrangements
    productive = dat.loc[dat['productive'] == 'T']
    counts = coll.Counter()
    reads = coll.defaultdict(fxn.nest_counter)
    # go through and get unique TCRs, with the most common TCR sequence to accompany them
    for row in productive.index:
        row_dat = productive.loc[row]
        vjcdr3 = row_dat['v_call'] + '|' + row_dat['j_call'] + '|' + row_dat['junction_aa']
        counts[vjcdr3] += 1
        reads[vjcdr3][row_dat['inferred_full_nt']] += 1
    #
    # Then go through those unique TCRs and output those with an above-threshold frequency to the outfile
    with opener('coll-' + input_args['in_file'], 'w') as out_file:
        out_file.write('\t'.join(out_headers) + '\n')
        line_count = 1
        for tcr in counts.most_common():
            # Cutoff if less frequent than the defined threshold
            if tcr[1] > input_args['threshold']:
                tcr_bits = tcr[0].split('|')
                if len(tcr_bits[2]) > 30:  # chuck out suspiciously long CDR3s missed by autoDCR
                    continue
                out_tcr = coll.defaultdict(str, {
                           'sequence_id': input_args['in_file'].split('.')[0] + '+' + str(line_count).zfill(5),
                           'v_call': tcr_bits[0],
                           'j_call': tcr_bits[1],
                           'junction_aa': tcr_bits[2],
                           'duplicate_count': str(tcr[1]),
                           'sequence': reads[tcr[0]].most_common(1)[0][0]
                           })
                out_file.write('\t'.join([out_tcr[x] for x in out_headers]) + '\n')
                line_count += 1



# LOOK FOR AN ATG AROUND 6-10 NT IN RUNNING UP TO THE START OF THE RC1 PRIMERS -- DEFINE THE EXACT TCR
# THEN CAN COLLAPSE EXACT SEQUENCES, MAP LENGTHS ETC

