from __future__ import print_function

import argparse
import sys
import khmer
import textwrap
import os

import pandas as pd
from khmer import Countgraph
from khmer.kfile import check_input_files
from khmer.khmer_args import (sanitize_help, KhmerArgumentParser)
from khmer._oxli.hashing import Kmer

from boink.analysis import find_N_most_abundant_kmers


def get_parser():
    epilog = """\

    """
    parser = KhmerArgumentParser(citations=['counting'])

    parser.add_argument('input_count_graph_filename', help='The name of the'
                        ' input k-mer countgraph file.')
    parser.add_argument('input_sequence_filenames', help='The name of the input'
                        ' FAST[AQ] sequence file.', nargs='+')
    parser.add_argument('-N', type=int, default=10000)
    parser.add_argument('-o', dest='output', type=argparse.FileType('w'), default=sys.stdout)
    
    return parser

def main():
    args = get_parser().parse_args()
    infiles = [args.input_count_graph_filename] + args.input_sequence_filenames
    for infile in infiles:
        check_input_files(infile, False)
    counts = khmer.load_countgraph(args.input_count_graph_filename)
    results = find_N_most_abundant_kmers(args.input_sequence_filenames,
                                         args.N, counts)

    results_df = pd.DataFrame({'kmer': [str(k) for k in results.keys()],
                               'count': [int(c) for c in results.values()]})
    results_df.sort_values(by='count', inplace=True, ascending=False)
    results_df.to_csv(args.output, index=False)

if __name__ == '__main__':
    main()
