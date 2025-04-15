# -*- coding: utf-8 -*-

"""
autoDCR

Automatic Decombinator: much like the regular DCR, this uses Aho-Corasick tries to locate short tag sequences
  in TCR sequencing data, and perform V/J annotation and subsequent CDR3 determination.

It's 'automatic' in the sense that it uses larger, unbiased tag set which tile across the entire V/J gene sets
  from IMGT. This allows comprehensive, allele-level resolution specification, at the cost of speed.

Looks for alpha and beta chain rearrangements simultaneously.

"""


import argparse
import gzip
import os
import collections as coll
from acora import AcoraBuilder
from time import time

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.2.7'
__author__ = 'Jamie Heather'


def args():
    """
    :return: argparse details of command line arguments passed to the script
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description='autoDCR v' + __version__ + ': find rearranged TCR sequences in HTS data, \n'
                    'using automatically generated tag sets.')
    # Add arguments
    parser.add_argument('-fq', '--fastq', type=str, required=True,
                        help='Correctly demultiplexed/processed FASTQ file containing TCR reads')

    parser.add_argument('-o', '--out_path', type=str, required=False,
                        help='Optionally specify path to output file.')

    parser.add_argument('-or', '--orientation', type=str, required=False, default="both",
                        help='Specify the orientation to search in (forward/reverse/both). Default = both')

    parser.add_argument('-sp', '--species', type=str, required=False, default="human",
                        help='Specify which species TCR repertoire the data consists of. Default = human')

    parser.add_argument('-dd', '--data_dir', type=str, required=False,
                        help="Optionally specify a path to a directory containing the required germline TCR data \n"
                             "(i.e. the 'X.fasta', 'X.tags', and 'X.translate' files)")

    parser.add_argument('-bc', '--barcoding', action='store_true', required=False,
                        help='Option to run autoDCR on barcoded libraries: a specified number of nucleotides will be '
                             'separated from the start of each read, for barcoding purposes.')

    parser.add_argument('-bl', '--bclength', type=int, required=False, default=42,
                        help='Length of barcode sequence, if applicable. Default is set to 42 bp.')

    parser.add_argument('-dl', '--deletion_limit', type=int, required=False, default=30,
                        help='Upper limit of allowable deletions for each of V/J in a recombination. Default = 30.')

    parser.add_argument('-cl', '--cdr3_limit', type=int, required=False, default=30,
                        help='Upper limit of allowable length of translated CDR3 junctions. Default = 30. '
                             'Set to 0 for no limit.')

    parser.add_argument('-jv', '--jump_values', action='store_true', required=False,
                        help="Optionally output the V and J 'jump' values, which record the position of the outermost "
                             "edges of the detected TCR rearrangements.")

    parser.add_argument('-ad', '--allele_discovery', action='store_true', required=False,
                        help='Run autoDCR in a mode that permits downstream attempts to discover novel TCR V/J alleles.'
                             '\nNote this is a highly experimental feature, and requires that tags have been generated '
                             'with default parameters (20 nt length, overlapping 10 nt). See README')

    parser.add_argument('-it', '--internal_translation', action='store_true', required=False,
                        help='Generate inferred full-length sequences off the most internal (CDR3-proximal) tags, '
                             'rather than the farthest (most distal) ones as per usual.'
                             '\nNote this is a highly experimental feature, but can be useful for applications where '
                             'users can reasonably expect variations indels/unexpected splicing in V/J genes.')

    parser.add_argument('-dt', '--dont_translate', action='store_true', required=False,
                        help='Stop the automatic translation of TCRs.')

    parser.add_argument('-dz', '--dont_gzip', action='store_true', required=False,
                        help='Stop the output FASTQ files automatically being compressed with gzip')

    return parser.parse_args()


def readfq(fastx_file):
    """
    readfq(file):Heng Li's Python implementation of his readfq function
    https://github.com/lh3/readfq/blob/master/readfq.py
    :param fastx_file: opened file containing fastq or fasta reads
    :yield: read id, read sequence, and (where available) read quality scores
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:
        if not last:  # the first record or a record following a fastq
            for l in fastx_file:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:], [], None  # This version takes the whole line (post '>')
        for l in fastx_file:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break

        else:  # this is a fastq record
            sequence, leng, seqs = ''.join(seqs), 0, []
            for l in fastx_file:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(sequence):  # have read enough quality
                    last = None
                    yield name, sequence, ''.join(seqs)  # yield a fastq record
                    break

            if last:  # reach EOF before reading enough quality
                yield name, sequence, None  # yield a fasta record instead
                break


def rev_comp(seq):
    """
    :param seq: Any DNA string, composed of just the four bases (upper or lower case)
    :return: The reverse complement of that DNA string (maintaining same case)

    """
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def get_gene(string):
    """
    :param string: a string containing an IMGT gene name - likely a FASTA header, or something derived from it
    :return: the IMGT gene/allele
    """
    return string.split('|')[1]


def build_trie(sequence_list):
    """
    :param sequence_list: A list of DNA sequences (tags) to search for
    :return: An Aho-Corasick trie, which can be used to search sequences downstream with all of the listed tags
    """

    trie_builder = AcoraBuilder()
    for sequence in sequence_list:
        trie_builder.add(sequence)

    return trie_builder.build()


def opener(file_name, open_mode):
    """
    :param file_name: Path to file to open
    :param open_mode: 'open' mode parameter, e.g. 'a', 'w', 'r'
    :return: The correct opener function (open or gzip.open)
    """
    if file_name.endswith('.gz'):
        return gzip.open(file_name, open_mode + 't')
    else:
        return open(file_name, open_mode)


def import_tcr_info(input_arguments):
    """
    Establishes global dictionaries which contain all the basic TCR information required for Decombining:
    - genes : All TCR V/J genes of the locus { gene_name: sequence }
    - tag_genes : The genes covered by that tag { tag_seq : gene_list }
    - tag_order : The position/index of each tag in the input tag file { tag_seq : index
    :param input_arguments: CLI arguments, containing a possible path to a data directoy
    :return: input_arguments (with potentially fixed data dir), plus all dictionaries are set up as global objects
    """

    globals()["genes"] = {}
    globals()["tag_genes"] = {}
    globals()["tag_order"] = {}

    d_dir = ''
    if 'data_dir' in input_arguments:
        if input_arguments['data_dir']:
            d_dir = input_arguments['data_dir']
            if not d_dir.endswith('/'):
                d_dir += '/'
                input_arguments['data_dir'] = d_dir

    # Get full gene sequences
    in_file_path = d_dir + input_arguments['species'] + '.fasta'
    if not os.path.exists(in_file_path):
        raise IOError("Cannot find input TCR FASTA file " + in_file_path + ".")

    with open(in_file_path, 'r') as in_file:
        for fasta_header, fasta_seq, quality in readfq(in_file):
            globals()["genes"][get_gene(fasta_header)] = fasta_seq.upper()

    # And also the tag information...
    in_file_path = d_dir + input_arguments['species'] + '.tags'
    if not os.path.exists(in_file_path):
        raise IOError("Cannot find input TCR TAG file " + in_file_path + ".")

    with open(in_file_path, 'r') as in_file:
        count = 0
        for line in in_file:
            tag_bits = line.rstrip().split('\t')
            globals()["tag_genes"][tag_bits[0]] = tag_bits[2]
            globals()["tag_order"][tag_bits[0]] = count

            count += 1

    # ... and built into Aho-Corasick tries, which will be used to search the data
    globals()["trie"] = build_trie(list(globals()["tag_genes"].keys()))

    # Also establish the global field of whether the allele discovery mode is enabled
    global discover
    if input_arguments['allele_discovery']:
        discover = True
    else:
        discover = False

    global internal
    if input_arguments['internal_translation']:
        internal = True
    else:
        internal = False

    # TODO add option to read in a pre-made trie?
    return input_arguments


def import_translate_info(input_arguments):
    """
    Establishes global dictionaries which contain all the information required for translating TCRs:
    - trans_pos : Position of conserved CDR3 junction defining residues { gene_name: index }
    - trans_res : Identity of conserved CDR3 junction defining residues { gene_name: residue_character }
    :param input_arguments: CLI input arguments
    :return: Nothing, all dictionaries are set up as global objects
    """

    globals()["trans_pos"] = {}
    globals()["trans_res"] = {}

    d_dir = ''
    if 'data_dir' in input_arguments:
        if input_arguments['data_dir']:
            d_dir = input_arguments['data_dir']

    # Get data
    in_file_path = d_dir + input_arguments['species'] + '.translate'
    if not os.path.exists(in_file_path):
        raise IOError("Cannot find input TCR translation file " + in_file_path + ".")

    with open(in_file_path, 'r') as in_file:
        for line in in_file:
            bits = line.rstrip().split('\t')
            globals()["trans_pos"][bits[0]] = int(bits[1])
            globals()["trans_res"][bits[0]] = bits[2]


def get_deletions(results_dict, input_arguments):
    """
    :param results_dict: The dictionary containing all of the TCR details discovered so far
    :param input_arguments: The argparse input arguments
    :return: the same results_dict, with additional information relating to the deletions discovered in the germline V/J
    """

    # TODO filter out cases where end of gene is beyond the limits of the read

    germlines = coll.defaultdict()
    for gene in ['V', 'J']:

        # Take the equivalent un-deleted sequence of that germline gene to compare against the recombined sequence
        full_germline = genes[results_dict[gene.lower() + '_call'].split(',')[0]]
        if gene == 'V':
            try:
                tag_start = full_germline.index(results_dict['inter_tag_seq'][:20])
                germlines[gene] = full_germline[tag_start:]
                recombined = results_dict['inter_tag_seq'][:len(germlines[gene])]
            except Exception:
                results_dict[gene.lower() + '_deletion_found'] = False
                break
        elif gene == 'J':
            try:
                tag_site = full_germline.index(results_dict['inter_tag_seq'][-20:])
                germlines[gene] = full_germline[:tag_site + 20]
                recombined = results_dict['inter_tag_seq'][-len(germlines[gene]):]
            except Exception:
                results_dict[gene.lower() + '_deletion_found'] = False
                break

        # Need to start at the right place and take sliding windows in the correction direction
        # V genes = Start at 3' and move 5' / J genes = start at 5' and move 3'
        position = {'V': -10, 'J': 0}
        increment = {'V': -1, 'J': 1}
        matched = False
        deletions = 0

        # Starting at the end of the gene, move the sliding window 1 nt away each iteration
        while not matched and deletions < len(recombined):
            if recombined[position[gene]:][:10] == germlines[gene][position[gene]:][:10]:
                matched = True
            else:
                position[gene] += increment[gene]
                deletions += 1

        if matched or deletions < input_arguments['deletion_limit']:
            results_dict[gene.lower() + '_deletions'] = deletions
            results_dict[gene.lower() + '_deletion_found'] = True
        else:
            results_dict[gene.lower() + '_deletion_found'] = False

    # Use the values determined above to pull out the 'insert' sequence, i.e. the non-template nt between ends of V/J
    if results_dict['v_deletion_found'] and results_dict['j_deletion_found']:
        its = results_dict['inter_tag_seq']
        try:
            # Use the length of the V (minus deletions) up to the 10-mer of post-deleted J sequence found above
            results_dict['insertion'] = its[len(germlines['V']) - results_dict['v_deletions']:its.index(
                germlines['J'][results_dict['j_deletions']:results_dict['j_deletions'] + 10])]
        except Exception:
            results_dict[gene.lower() + '_deletion_found'] = False

    return results_dict


def dcr(sequence, quality, passed_input_arguments):
    """
    Core wrapper function that performs the actual decombining
    :param sequence: DNA string
    :param quality: Corresponding quality scores (if available - padded with spaces if not)
    :param passed_input_arguments: passing the input arguments dict along for get_deletions
    :return: results dictionary containing the details of any discovered recombined TCR
    """

    # Search the read using the Aho-Corasick trie, use this to generate an inter-tag sequence
    check = find_tag_hits(sequence)

    if check:
        check['inter_tag_seq'] = sequence[check['v_tag_position']:check['j_tag_position'] + 20]
        if quality:
            check['inter_tag_qual'] = quality[check['v_tag_position']:check['j_tag_position'] + 20]

        get_deletions(check, passed_input_arguments)

    return check


def find_tag_hits(sequence):
    """
    Does the actual tag searching
    :param sequence: DNA seq to search for tag hits (i.e. sequences matching tags corresponding to TCR genes)
    :return: defaultdict containing fields corresponding to the most distal/most unique V/J tag hits
    """

    # Search the read using the relevant Aho-Corasick trie, and pull out the genes those tags come from
    # call, tag_index, match, jump, i, len_check = '', '', ['', ''], '', '', ''
    out_dict = coll.defaultdict()
    check = globals()['trie'].findall(sequence)

    found = 0
    for gene_type in ['V', 'J']:
        specific_check = [x for x in check if gene_type in tag_genes[x[0]].split('*')[0]]

        if specific_check:
            found += 1
            gene_lists = [tag_genes[x[0]] for x in specific_check]

            # Go through and find the minimal set of genes that are featured in the maximum number of tags
            gene_counts = coll.Counter()
            for gene_list in [x.split(',') for x in list(gene_lists)]:
                for g in gene_list:
                    gene_counts[g] += 1

            for i in range(len(specific_check), 0, -1):
                hits = [x for x in gene_counts if gene_counts[x] == i]
                if len(hits) > 0:
                    hits.sort()
                    out_dict[gene_type.lower() + '_call'] = ','.join(hits)
                    break

            out_dict[gene_type.lower() + '_number_tags_matches'] = i
            out_dict[gene_type.lower() + '_number_tags_total'] = len(specific_check)

            # Then go through the matched tag and find the relevant hit that uses the tag combination found
            if gene_type == 'V':
                if not internal:
                    tag_order = [x for x in specific_check]
                else:
                    tag_order = [x for x in specific_check][::-1]
            elif gene_type == 'J':
                if not internal:
                    tag_order = [x for x in specific_check][::-1]
                else:
                    tag_order = [x for x in specific_check]
            else:
                raise ValueError("Unexpected gene type during translation. ")

            # Some genes are uniquely identifiable only through their unique appearance in the intersection of many tags
            # Need to find the most-distal tag that covers that gene/all of those genes
            # Use reference positions from the only OR first listed gene
            call_bits = out_dict[gene_type.lower() + '_call'].split(',')

            for match in tag_order:
                tag_genes_sep = tag_genes[match[0]].split(',')

                if all(g in tag_genes_sep for g in call_bits):

                    out_dict[gene_type.lower() + '_jump'] = genes[hits[0]].index(match[0])

                    # Need to add on the J tag length to make the inter-tag sequence include the J tag
                    if gene_type == 'J':
                        out_dict[gene_type.lower() + '_jump'] += len(match[0])

                    break  # Only need to find the outermost jump values

                else:
                    # Of the V|J tags found, there are no alleles represented among all and thus no unambig call
                    pass
                    # TODO if desired not-rearranged-but-TCR-containing reads could be output here

            # Check whether allele discovery mode has been established
            if 'discover' in globals():
                if discover:
                    # Only passrd through tags that correspond to the called gene
                    gene_specific_tag_order = [x for x in tag_order if out_dict[gene_type.lower() + '_call']
                                               in tag_genes[x[0]].split(',')]

                    # Can only perform allele discovery sampling if there's a gene call/remaining tags
                    if len(gene_specific_tag_order) > 0:
                        out_dict[gene_type.lower() + '_mismatches'] = noncontiguous_tag_check(
                            sequence, gene_specific_tag_order, 10, out_dict[gene_type.lower() + '_call'])

            out_dict[gene_type.lower() + '_tag_seq'] = match[0]
            out_dict[gene_type.lower() + '_tag_position'] = match[1]

    # If both a V and a J gene have been called, return the output
    if found == 2:
        return out_dict


def noncontiguous_tag_check(sequence, tag_hits, win_len, gene):
    """
    Used in trying to identify potential novel alleles, based on a break in consecutive tag hits
    :param sequence: read being decombined
    :param tag_hits: the list of V/J tag hits produced by find_tag_hits
    :param win_len: length of expected window between adjacent tags
    :param gene: gene used
    :return: dict detailing the detected noncontiguous mismatch, i.e. the tags before/after and sequence in between
    """

    # Account for fact J gene tag orders are running in descending order
    if gene[3] == 'J':
        tag_hits = tag_hits[::-1]

    # Then simply find 'missing' tags, based on what tag positions we'd expect to see, given a certain tag window length
    tag_positions = [x[1] for x in tag_hits]

    non_contig = [x for x in range(min(tag_positions), max(tag_positions), win_len) if x not in tag_positions]

    # TODO add check here that all tags map to the same gene?

    if not non_contig:
        return ''
    else:
        # Find the substring that's not covered by the tags, throwing it out if the tags before/after aren't in sync
        try:
            # Find the tags and positions in the read
            tag1_seq, tag1_index = [x for x in tag_hits if x[1] == non_contig[0] - win_len][0]
            tag2_seq, tag2_index = [x for x in tag_hits if x[1] == non_contig[-1] + win_len][0]

            # And the corresponding parts in the germline gene
            tag1_gl_index = genes[gene.split(',')[0]].index(tag1_seq)
            tag2_gl_index = genes[gene.split(',')[0]].index(tag2_seq)

        except Exception:
            return ''

        # Otherwise report back all the relevant bracketing and intervening sequence/tag information
        mismatch_dict = {
            'tag1_index': tag1_index,
            'tag2_index': tag2_index,
            'tag1_seq': tag1_seq,
            'tag2_seq': tag2_seq,
            'tag1_gl_index': tag1_gl_index,
            'tag2_gl_index': tag2_gl_index,
            'intervening': sequence[tag1_index + len(tag1_seq):tag2_index]
        }
        return mismatch_dict


def translate(nt_seq):
    """
    :param nt_seq: Nucleotide sequence to be translated
    :return: corresponding amino acid sequence
    """

    aa_seq = ''
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3]
        if len(codon) == 3:
            try:
                aa_seq += codons[codon]
            except Exception:
                raise IOError("Cannot translate codon: " + codon)

    return aa_seq


codons = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
          'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
          'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
          'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
          'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
          'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
          'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
          'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
          'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
          'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
          'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
          'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
          'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
          'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
          'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
          'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}


def sort_read_bits(whole_read, whole_qual, input_arguments):
    """
    :param whole_read: The entire read, from start to finish
    :param whole_qual: The quality scores of that entire reads
    :param input_arguments: The argparse input arguments, to determine whether/how to chop up the read
    :return: 4 str: the TCR-containing portion of the read, it's quality, and the barcode-containing read/its quality
    """
    if input_arguments['barcoding']:
        tcr_read = whole_read[input_arguments['bclength']:].upper()
        tcr_qual = whole_qual[input_arguments['bclength']:]
        bc_read = whole_read[:input_arguments['bclength']].upper()
        bc_qual = whole_qual[:input_arguments['bclength']]
        return tcr_read, tcr_qual, bc_read, bc_qual
    else:
        return whole_read.upper(), whole_qual, '', ''


def tcr_search(tcr_read, tcr_qual, input_arguments, headers):
    """
    Tries to look for TCRs in all of the requested frames
    :param tcr_read: Str containing the portion of the read that should/might contain rearrangements
    :param tcr_qual: Str containing corresponding quality scores
    :param input_arguments: Dict of argparse input arguments, containing the needed frame
    :param headers: List of fields to be included in final output, even if only as empty fields
    :return: Dict of rearranged TCR properties, if it finds one
    """

    search = ''
    if input_arguments['orientation'] == 'forward' or input_arguments['orientation'].lower() == 'f':
        search = dcr(tcr_read.upper(), tcr_qual, input_arguments)
        if search:
            search['rev_comp'] = 'F'

    elif input_arguments['orientation'] == 'reverse' or input_arguments['orientation'].lower() == 'r':
        search = dcr(rev_comp(tcr_read.upper()), tcr_qual[::-1], input_arguments)
        if search:
            search['rev_comp'] = 'T'

    elif input_arguments['orientation'] == 'both' or input_arguments['orientation'].lower() == 'b':
        # If searching in both orientations, need to open to the possibility of finding TCRs in both
        search_f = dcr(tcr_read.upper(), tcr_qual, input_arguments)
        search_r = dcr(rev_comp(tcr_read.upper()), tcr_qual[::-1], input_arguments)

        if search_f and not search_r:
            search = search_f
            search['rev_comp'] = 'F'
        elif not search_f and search_r:
            search = search_r
            search['rev_comp'] = 'T'
        # TODO potentially could add 'elif search_f and search_r' here: unlikely to have data with bidirectional TCRs?

    else:
        raise IOError("Incorrect orientation argument provided: " + input_arguments['orientation'] + "\n"
                      "Please specify one of the three options: forward/reverse/both - or f/r/b.")

    if search:
        search['sequence'] = tcr_read

        if not input_arguments['dont_translate']:
            search = find_cdr3(search)

        # Finally pad dict with empty values for all columns required in final output document (AIRR community format)
        for field in headers:
            if field not in search:
                search[field] = ''

    return search


def find_cdr3(tcr):
    """
    :param tcr: Dict of recombination properties, as produced by the dcr and tcr_search functions
    :return: Same dict, with corresponding translated/CDR3 search properties entered in
    """
    # Need to check for presence of both jump values
    # Can lack one whilst still having a detected rearrangement if say multiple Vs/Js in one read
    if 'v_jump' in tcr and 'j_jump' in tcr:
        tcr['inferred_full_nt'] = genes[tcr['v_call'].split(',')[0]][:tcr['v_jump']] \
                                  + tcr['inter_tag_seq'] + genes[tcr['j_call'].split(',')[0]][tcr['j_jump']:]
        tcr['inferred_full_aa'] = translate(tcr['inferred_full_nt'])

        # If multiple genes detected, just arbitrarily pick one to use for translation parameters
        if ',' in tcr['v_call']:
            translate_v = tcr['v_call'].split(',')[0]
        else:
            translate_v = tcr['v_call']

        if ',' in tcr['j_call']:
            translate_j = tcr['j_call'].split(',')[0]
        else:
            translate_j = tcr['j_call']

        tcr['junction_aa'] = tcr['inferred_full_aa'][trans_pos[translate_v]:trans_pos[translate_j] + 1]

        # Assume productivity, and remove it as applicable, checking for the various required parameters
        # Note that it's possible for multiple different reasons to be non-productive to be true, in different combos
        tcr['productive'] = 'T'
        tcr['vj_in_frame'] = 'T'

        # Check whether it's feasibly in frame
        if (len(tcr['inferred_full_nt']) - 1) % 3 != 0:
            tcr['productive'] = 'F'
            tcr['vj_in_frame'] = 'F'

        # Need to account for cases where there is no detectable/valid CDR3
        if tcr['junction_aa']:
            # Check for stop codons...
            if '*' in tcr['inferred_full_aa']:
                tcr['stop_codon'] = 'T'
                tcr['productive'] = 'F'
            else:
                tcr['stop_codon'] = 'F'

            # ... and the conserved V gene residue at the right position...
            if tcr['junction_aa'][0] != trans_res[translate_v]:
                tcr['conserved_c'] = 'F'
                tcr['productive'] = 'F'
            else:
                tcr['conserved_c'] = 'T'

            # ... and same for the J gene...
            if tcr['junction_aa'][-1] != trans_res[translate_j]:
                tcr['conserved_f'] = 'F'
                tcr['productive'] = 'F'
            else:
                tcr['conserved_f'] = 'T'

            # And check whether the CDR3 falls within the expected length range
            if input_args['cdr3_limit'] > 0:
                if len(tcr['junction_aa']) <= input_args['cdr3_limit']:
                    tcr['cdr3_in_limit'] = 'T'
                else:
                    tcr['cdr3_in_limit'] = 'F'
                    tcr['junction_aa'] = ''
                    tcr['productive'] = 'F'

            else:
                tcr['cdr3_in_limit'] = ''


        else:
            tcr['productive'] = 'F'

        # Cryptic splices in leader processing can result in a variety of translation issues, let's catch those
        if tcr['inferred_full_aa'] and not tcr['junction_aa'] and tcr['inter_tag_seq']:
            tcr = attempt_salvage_irregular_cdr3s(tcr, translate_v, translate_j)

        if tcr['productive'] == 'F':

            if tcr['junction_aa']:

                # Additional check to try to salvage recovery of CDR3s in frame-shifted receptors
                # Works off scanning from conserved F (if present) and scanning N-wards for the next C
                if tcr['stop_codon'] == 'T' and tcr['conserved_f'] == 'F':
                    tcr = attempt_salvage_irregular_cdr3s(tcr, translate_v, translate_j)

                if 'non_productive_junction_aa' not in tcr:
                    tcr['non_productive_junction_aa'] = tcr['junction_aa']
                tcr['junction_aa'] = ''

        else:
            tcr['junction'] = tcr['inferred_full_nt'][
                              trans_pos[translate_v] * 3:(trans_pos[translate_j] * 3) + 2]

    return tcr


def attempt_salvage_irregular_cdr3s(tcr_dict, v_translate, j_translate):
    """
    Try to find CDR3s in irregular rearrangements (particularly necessary when using irregular V-REGION references)
    :param tcr_dict: Dictionary describing rearrangement under consideration
    :param v_translate: str ID of V gene being used for translation purposes (maybe differ if >1 detected)
    :param j_translate: str ID of V gene being used for translation purposes (maybe differ if >1 detected)
    :return: tcr_dict, hopefully updated with any CDR3s the regular assumptions failed to catch now annotated
    """
    potential_frames = []
    potential_translations = []

    for frame in range(3):
        potential_translation = translate(tcr_dict['inferred_full_nt'][frame:])

        if potential_translation[trans_pos[j_translate]] == trans_res[j_translate]:
            potential_frames.append(frame)
            potential_translations.append(potential_translation)

    # If there's two frames, one which has stops and one which doesn't, keep the one without
    if len(potential_frames) == 2:
        discard = []
        for i in [0, 1]:
            if '*' in potential_translations[i]:
                discard.append(i)
        if len(discard) == 1:
            potential_frames.pop(discard[0])

    # If there's a single frame in which there's a conserved F residue at the appropriate location...
    if len(potential_frames) == 1:
        # ... look upstream until it finds a C (within a typical distance)
        potential_translation = translate(tcr_dict['inferred_full_nt'][potential_frames[0]:])
        j_f_pos = trans_pos[j_translate]
        for i in range(8, 21):
            potential_junction_aa = potential_translation[j_f_pos - i:j_f_pos + 1]
            if potential_junction_aa[0] == trans_res[v_translate]:
                tcr_dict['non_productive_junction_aa'] = potential_junction_aa
                tcr_dict['inferred_full_aa'] += ',' + translate(
                    tcr_dict['inferred_full_nt'][potential_frames[0]:])
                tcr_dict['conserved_f'] = 'T'
                break

    return tcr_dict


out_headers = ['sequence_id', 'v_call', 'd_call', 'j_call', 'junction_aa', 'duplicate_count', 'sequence',
               'junction', 'rev_comp', 'productive', 'sequence_aa',
               'inferred_full_nt', 'inferred_full_aa', 'non_productive_junction_aa',
               'vj_in_frame', 'stop_codon', 'conserved_c', 'conserved_f', 'cdr3_in_limit',
               'inter_tag_seq', 'inter_tag_qual', 'umi_seq', 'umi_qual',
               'sequence_alignment', 'germline_alignment', 'v_cigar', 'd_cigar', 'j_cigar']


if __name__ == '__main__':

    # Determine the requested input/output parameters
    input_args = vars(args())
    input_args = import_tcr_info(input_args)

    if not input_args['dont_translate']:
        import_translate_info(input_args)

    if discover:
        out_headers += ['v_mismatches', 'j_mismatches']
        input_args['jump_values'] = True

    if input_args['jump_values']:
        out_headers += ['v_jump', 'j_jump']

    counts = coll.Counter()
    start = time()

    # TODO sanity/presence check the input FQ (including a length check - give warning if too short)
    # Determine where to save the results
    analysis_name = input_args['fastq'].split('/')[-1].split('.')[0]
    if discover:
        analysis_name += '_infer-alleles'
    out_file_name = analysis_name + '.tsv'

    if 'out_path' in input_args:
        if input_args['out_path']:
            if input_args['out_path'].endswith('/'):
                out_file_name = input_args['out_path'] + analysis_name + '.tsv'
            else:
                out_file_name = input_args['out_path'].replace('.tsv.', '').replace('.gz', '') + '.tsv'

    if not input_args['dont_gzip']:
        out_file_name += '.gz'

    with opener(input_args['fastq'], 'r') as in_file, opener(out_file_name, 'w') as out_file:

        # Initialise the output file with the header
        out_file.write('\t'.join(out_headers) + '\n')
        out_str = []
        n_skip = 0
        for read_id, seq, qual in readfq(in_file):

            # Pad empty quality scores for FASTA files
            if not qual:
                qual = ' ' * len(seq)

            counts['reads'] += 1

            # Figure out the relevant parts of the read for decombining, then search
            read, real_qual, bc, bc_qual = sort_read_bits(seq, qual, input_args)
            if 'N' in read:
                n_skip += 1
                continue
            tcr_check = tcr_search(read, real_qual, input_args, out_headers)
            if tcr_check:

                if input_args['barcoding']:
                    tcr_check['umi_seq'] = bc
                    tcr_check['umi_qual'] = bc_qual

                counts['rearrangements'] += 1
                tcr_check['sequence_id'] = read_id

                line_out = '\t'.join([str(tcr_check[x]) for x in out_headers])

                if discover:
                    if tcr_check['v_mismatches'] or tcr_check['j_mismatches']:
                        counts['mismatched_germlines'] += 1

                out_str.append(line_out)

                # Bulk write the results out once there's a sufficient chunk (to prevent this getting too big in memory)
                if len(out_str) % 10000 == 0:
                    out_file.write('\n'.join(out_str) + '\n')
                    out_str = []

        # Then write out any leftover calls
        out_file.write('\n'.join(out_str))
        # TODO fix duplicate counts (if desired?)

    end = time()
    time_taken = end - start
    print("Took", str(round(time_taken, 2)), "seconds")
    print("Found", str(counts['rearrangements']), "rearranged TCRs in", str(counts['reads']), "reads")
    if discover:
        print("Of these,", str(counts['mismatched_germlines']), "showed discontinuous tag matches "
              "and were kept aside for inference of potential new alleles.")
    print("\t" + str(n_skip) + " reads skipped due to containing N calls?")
    # TODO sort summary output (maybe into YAML?)
