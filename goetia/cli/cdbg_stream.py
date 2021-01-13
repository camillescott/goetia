#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : cdbg_stream.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 11.03.2020

from goetia import libgoetia

from goetia.cdbg import (compute_connected_component_callback,
                         compute_unitig_fragmentation_callback,
                         write_cdbg_metrics_callback,
                         write_cdbg_callback)
from goetia.dbg import get_graph_args, process_graph_args
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.processors import AsyncSequenceProcessor
from goetia.messages import (Interval, SampleStarted, SampleFinished, Error, AllMessages)
from goetia.metadata import CUR_TIME
from goetia.serialization import cDBGSerialization

from goetia.cli.args import get_output_interval_args, print_interval_settings
from goetia.cli.runner import CommandRunner

import curio

import os
import sys


class cDBGRunner(CommandRunner):

    def __init__(self, parser):
        get_graph_args(parser)
        get_cdbg_args(parser)
        get_output_interval_args(parser)

        group = get_fastx_args(parser)
        group.add_argument('-o', dest='output_filename', default='/dev/stdout')
        group.add_argument('-i', '--inputs', dest='inputs', nargs='+', required=True)

        parser.add_argument('--echo', default=None,
                            help='echo all events to the given file.')
        parser.add_argument('--curio-monitor', default=False, action='store_true',
                            help='Run curio kernel monitor for async debugging.')

        super().__init__(parser)

    def postprocess_args(self, args):
        process_graph_args(args)
        process_cdbg_args(args)

    def setup(self, args):
        os.makedirs(args.results_dir, exist_ok=True)

        self.dbg_t       = args.graph_t
        self.hasher      = args.hasher_t(args.ksize)
        self.storage     = args.storage.build(*args.storage_args)
        self.dbg         = args.graph_t.build(self.storage, self.hasher)

        self.cdbg_t      = libgoetia.cdbg.cDBG[type(self.dbg)]

        self.compactor_t = libgoetia.cdbg.StreamingCompactor[type(self.dbg)]

        self.compactor = self.compactor_t.Compactor.build(self.dbg)

        if args.normalize:
            self.file_processor = self.compactor_t.NormalizingCompactor[FastxReader].build(self.compactor,
                                                                                           args.normalize,
                                                                                           args.fine_interval,
                                                                                           args.medium_interval,
                                                                                           args.coarse_interval)
        else:
            self.file_processor = self.compactor_t.Processor.build(self.compactor,
                                                                   args.fine_interval,
                                                                   args.medium_interval,
                                                                   args.coarse_interval)
        
        # Iterator over samples (pairs or singles, depending on pairing-mode)
        sample_iter = iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names)
        # AsyncSequenceProcessor does event management and callback for the FileProcessors
        self.processor = AsyncSequenceProcessor(self.file_processor, sample_iter, args.echo)
        # Subscribe a listener to the FileProcessor producer
        self.worker_listener = self.processor.add_listener('worker_q', 'cdbg.consumer')

        #
        # Register callbacks for data outputs.
        # Track a list of files that need to be closed with a ]
        # when we're done.
        # 
        self.to_close = []

        if args.track_cdbg_stats:
            self.worker_listener.on_message(Interval,
                                            write_cdbg_metrics_callback,
                                            self.compactor,
                                            args.track_cdbg_stats)
            self.to_close.append(args.track_cdbg_stats)


        if args.track_cdbg_unitig_bp:
            if args.unitig_bp_bins is None:
                bins = [args.ksize, 100, 200, 500, 1000]
            else:
                bins = args.unitig_bp_bins
            
            self.worker_listener.on_message(Interval,
                                            compute_unitig_fragmentation_callback,
                                            self.cdbg_t,
                                            self.compactor.cdbg,
                                            args.track_cdbg_unitig_bp,
                                            bins)
            self.to_close.append(args.track_cdbg_unitig_bp)


        if args.track_cdbg_components:
            self.worker_listener.on_message(Interval,
                                            compute_connected_component_callback,
                                            self.cdbg_t,
                                            self.compactor.cdbg,
                                            args.track_cdbg_components,
                                            args.component_sample_size)
            self.to_close.append(args.track_cdbg_components)

        if args.save_cdbg:
            for cdbg_format in args.save_cdbg_format:
                self.worker_listener.on_message(Interval,
                                                write_cdbg_callback,
                                                args.save_cdbg,
                                                cdbg_format)
                self.worker_listener.on_message(SampleFinished,
                                                write_cdbg_callback,
                                                args.save_cdbg,
                                                cdbg_format)

        # Close all files when done
        async def close_files(msg, files):
            for file_name in files:
                async with curio.aopen(file_name, 'a') as fp:
                    await fp.write('\n]\n')

        self.worker_listener.on_message(SampleFinished, close_files, self.to_close)

        #
        # Regular diagnostics output
        # 

        def info_output(msg):
            info = f'{msg.msg_type}: {getattr(msg, "state", "")}'\
                   f'\n\tSample:    {msg.sample_name}'\
                   f'\n\tSequences: {msg.t}'
            if msg.msg_type == 'Error':
                info += f'\n\tError: {msg.error}'

            print(info, file=sys.stderr)

        self.worker_listener.on_message(AllMessages, info_output)
                           
    def execute(self, args):
        curio.run(self.processor.start, with_monitor=args.curio_monitor)

    def teardown(self):
        pass


def get_cdbg_args(parser):
    default_prefix = 'goetia.build-cdbg.' + CUR_TIME
    parser.default_prefix = default_prefix
    group = parser.add_argument_group('cDBG')

    group.add_argument('--results-dir',
                        default=default_prefix)

    group.add_argument('--normalize',
                        type=int,
                        nargs='?',
                        const=10)

    group.add_argument('--save-cdbg',
                        metavar='PREFIX.<format>',
                        nargs='?',
                        const='goetia.cdbg.graph',
                        help='Save a copy of the cDBG.')
    group.add_argument('--save-cdbg-format',
                        nargs='+',
                        choices=cDBGSerialization.FORMATS,
                        default=['gfa1'])
    group.add_argument('--cdbg-tick',
                       type=int,
                       default=10,
                       help='Save every N interval ticks.')

    group.add_argument('--track-cdbg-metrics',
                        metavar='FILE_NAME.json',
                        nargs='?',
                        const='goetia.cdbg.stats.json',
                        help='Output basic cDBG metrics.')
    group.add_argument('--cdbg-metrics-tick',
                       type=int,
                       default=5,
                       help='Output every N interval ticks.')

    group.add_argument('--track-cdbg-components',
                        metavar='FILE_NAME.json',
                        nargs='?',
                        const='goetia.cdbg.components.json',
                        help='Save the distribution of component sizes.')
    group.add_argument('--component-sample-size',
                        type=int,
                        default=10000,
                        help='Number of components to sample for size.')
    group.add_argument('--cdbg-components-tick',
                       type=int,
                       default=5,
                       help='Sample and save distribution every N interval ticks.')

    group.add_argument('--track-unitig-bp',
                        metavar='FILENAME.json',
                        nargs='?',
                        const='goetia.cdbg.unitigs.bp.json',
                        help='Track the distribution of unitig sizes.')
    group.add_argument('--unitig-bp-bins',
                        nargs='+',
                        type=int,
                        help='Bin sizes of distribution.')
    group.add_argument('--unitig-bp-tick',
                       type=int,
                       default=10)

    group.add_argument('--validate',
                        metavar='FILENAME.csv',
                        nargs='?',
                        const='goetia.cdbg.validation.csv')

    return group


def process_cdbg_args(args):

    def join(p):
        return p if p is None else os.path.join(args.results_dir, p)

    args.track_cdbg_stats =      join(args.track_cdbg_stats)
    args.track_cdbg_components = join(args.track_cdbg_components)
    args.save_cdbg =             join(args.save_cdbg)
    args.track_cdbg_unitig_bp =  join(args.track_cdbg_unitig_bp)


def print_cdbg_args(args):
    print('* cDBG Params', file=sys.stderr)
    print('* Directory: ', args.results_dir, file=sys.stderr)
    if args.save_cdbg:
        print('* Saving cDBG every {0} sequences with file prefix {1}'.format(args.coarse_interval,
                                                                              args.save_cdbg),
              file=sys.stderr)
        print('* cDBG save formats: {0}'.format(', '.join(args.save_cdbg_format)))
    if args.track_cdbg_stats:
        print('* Tracking cDBG stats and reporting every {0} sequences'.format(args.fine_interval),
              file=sys.stderr)
        print('* Saving tracking information to', args.track_cdbg_stats, file=sys.stderr)
    if args.track_cdbg_history:
        print('* Tracking cDBG history and saving to', args.track_cdbg_history, file=sys.stderr)
    if args.validate:
        print('* cDBG will be validated on completion and results saved to', args.validate,
              file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)


