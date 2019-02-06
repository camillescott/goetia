#!/usr/bin/env python

import argparse
import sys

from boink.args import build_dBG_args, add_output_interval_args
from boink.dbg import make_dBG
from boink.compactor import make_streaming_compactor
from boink.processors import make_decision_node_processor

def parse_args():
    parser = build_dBG_args()
    add_output_interval_args(parser)
    parser.add_argument('-o', dest='output_filename', default='/dev/stdout')
    parser.add_argument('-i', dest='inputs', nargs='+', default=['/dev/stdin'])
    parser.add_argument('--output-interval', type=int, default=10000)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    graph = make_dBG(args.ksize,
                     args.max_tablesize,
                     args.n_tables,
                     storage='_' + args.storage_type)
    compactor = make_streaming_compactor(graph)
    processor = make_decision_node_processor(compactor,
                                             args.output_filename,
                                             args.fine_interval,
                                             args.medium_interval,
                                             args.coarse_interval)
    for filename in args.inputs:
        processor.process(filename)

    if args.savegraph:
        graph.save(args.savegraph)


if __name__ == '__main__':
    main()
