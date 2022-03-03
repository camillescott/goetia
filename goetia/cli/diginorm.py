#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : diginorm.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2021

import sys
import time

from goetia.filters import DiginormFilter
from goetia.dbg import get_graph_args, process_graph_args
from goetia.cli.args import get_output_interval_args
from goetia.cli.runner import CommandRunner
from goetia.cli.status import StatusOutput
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.processors import JSONStreamWriter
from goetia.storage import get_storage_args, process_storage_args
from goetia.utils import Counter, time_iterable

import blessings


desc = '''
{term.italic}filter module: {term.normal}diginorm filter

    Filter the input reads with digital (in silico) normalization
    at the specified median count cutoff.
'''

class DiginormFilterRunner(CommandRunner):

    def __init__(self, parser):
        get_storage_args(parser, default='NibbleStorage')
        get_graph_args(parser)
        get_output_interval_args(parser)

        group = get_fastx_args(parser)
        group.add_argument('-o', dest='output_filename', default='/dev/stdout')
        group.add_argument('-i', '--inputs', dest='inputs', nargs='+', required=True)
        group.add_argument('--quiet', action='store_true', default=False)

        parser.add_argument('-C', '--cutoff', type=int, default=10)
        parser.add_argument('--metrics', nargs='?', const='goetia.diginorm.metrics.json')
        
        super().__init__(parser, description=desc)

    def postprocess_args(self, args):
        process_graph_args(args)

    def setup(self, args):
        if args.metrics:
            self.metrics_stream = JSONStreamWriter(args.metrics)
        self.dbg_t       = args.graph_t
        self.hasher      = args.hasher_t(args.ksize)
        self.storage     = args.storage.build(*args.storage_args)
        self.dbg         = args.graph_t.build(self.storage, self.hasher)
        self.filter_t    = DiginormFilter[self.dbg_t]
        self.filter      = self.filter_t.Filter.build(self.dbg,
                                                      args.cutoff)

        self.processor = self.filter_t.Processor.build(self.filter.__smartptr__(),
                                                       args.output_filename,
                                                       args.interval)
        self.status = StatusOutput(['unique k-mers', 'sequences passed', 'estimated fp-rate'], term=self.term)

    def execute(self, args):
        warned_fp = False
        for (sample, name), sample_start_time, _ in \
            time_iterable(iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names)):

            self.status.start_sample(name, sample)
            for (n_seqs, stream_time, n_skipped), \
                 interval_start_time, interval_elapsed_time in \
                 time_iterable(self.processor.chunked_process(*sample)):
                
                n_passed = self.processor.n_passed()
                p_passed = n_passed / n_seqs
                passed = f'{n_passed:,} ({p_passed * 100.0:.1f}%)'
                fp_rate = self.dbg.estimated_fp()
                self.status.update(stream_time, n_seqs, interval_elapsed_time, self.dbg.n_unique(), passed, fp_rate)

                if args.metrics:
                    self.metrics_stream.write({'t': stream_time,
                                               'seq_t': n_seqs,
                                               'rt_elapsed_interval': interval_elapsed_time,
                                               'n_seqs_passed': n_passed,
                                               'p_seqs_passed': p_passed,
                                               'n_unique_kmers': self.dbg.n_unique(),
                                               'estimated_fp': fp_rate})
                
                if fp_rate >= 0.8 and warned_fp is False:
                    self.status.message(f'WARNING: false positive rate has surpassed 0.8 ({fp_rate:.7f})')
                    warned_fp = True


            sample_elapsed_time = time.perf_counter() - sample_start_time
            self.status.finish_sample(sample_elapsed_time)

    def teardown(self):
        pass
