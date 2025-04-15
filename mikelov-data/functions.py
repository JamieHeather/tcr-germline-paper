import pandas as pd
import os
import datetime
import numpy as np
import collections as coll


def fastqfy(header, sequence, quality):
    """
    """
    return "@" + header + "\n" + sequence + "\n+\n" + quality + '\n'


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
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fastx_file:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break

            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

def nest():
    """
    Create nested defaultdicts
    """
    return coll.defaultdict(list)


def nest_counter():
    """
    Create nested counters
    """
    return coll.Counter()


def list_to_df(list_o_lists, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param list_o_lists: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    df = pd.DataFrame(list_o_lists)
    df = df.rename(index=str, columns=dict(zip(range(len(headers)), headers)))
    df = df.sort_values(by=[headers[0]])
    if rename:
        df = df.set_index(headers[0], drop=True)
    return df


def plot_dir(plot_dir_suffix):
    """
    :return: The path to the plots directory subfolder for results plotted on this day (creating it if needed)
    """
    plot_dir = base_plot_dir + today() + '-' + plot_dir_suffix + '/'
    make_check_dir(plot_dir)
    return plot_dir


def make_check_dir(dir_path):
    """
    :param dir_path: A directory to check whether it exists; if not, make it
    """
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    return dir_path


def today():
    """
    :return: Today's day, in ISO format
    """
    return datetime.datetime.today().date().isoformat()


base_plot_dir = '../Plots/'


def rev_comp(seq):
    """
    :param seq: Any DNA string, composed of just the four bases (upper or lower case)
    :return: The reverse complement of that DNA string (maintaining same case)
    """
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]