#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : solid_filter.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 12.03.2020

import sys
import time

from goetia.filters import SolidFilter
from goetia.dbg import get_graph_args, process_graph_args
from goetia.cli.args import get_output_interval_args
from goetia.cli.runner import CommandRunner
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.storage import get_storage_args, process_storage_args
from goetia.utils import Counter

import blessings


desc = '''
{term.italic}filter module: {term.normal}solid filter

    Filter the input sequences by whether a specified proportion of their
    k-mers exceeds a specified minimum count threshold.
'''

class SolidFilterRunner(CommandRunner):

    def __init__(self, parser):
        get_storage_args(parser, default='NibbleStorage')
        get_graph_args(parser)
        get_output_interval_args(parser)
        group = get_fastx_args(parser)
        group.add_argument('-o', dest='output_filename', default='/dev/stdout')
        group.add_argument('-i', '--inputs', dest='inputs', nargs='+', required=True)
        parser.add_argument('-C', '--solid-min-count', type=int, default=1)
        parser.add_argument('-P', '--solid-min-proportion', type=float, default=.75)

        super().__init__(parser, description=desc)

    def postprocess_args(self, args):
        process_graph_args(args)

    def setup(self, args):
        self.dbg_t       = args.graph_t
        self.hasher      = args.hasher_t(args.ksize)
        self.storage     = args.storage.build(*args.storage_args)
        self.dbg         = args.graph_t.build(self.storage, self.hasher)
        self.filter_t    = SolidFilter[self.dbg_t]
        self.solid_filter = self.filter_t.Filter.build(self.dbg,
                                                       args.solid_min_proportion,
                                                       args.solid_min_count)

        self.processor = self.filter_t.Processor.build(self.solid_filter.__smartptr__(),
                                                       args.output_filename,
                                                       args.interval)
        self.status = SolidFilterRunner.StatusOutput(term=self.term)

    class StatusOutput:

        def __init__(self, term=None, file=sys.stderr):
            self.file = file
            self.term = term or blessings.Terminal()

        def print(self, *args, **kwargs):
            print(*args, **kwargs, file=self.file)

        def start_sample(self, sample_name, file_names):
            self.print(f'{self.term.italic}Begin sample: {self.term.normal}{sample_name}')
            files = '\n'.join(['    ' + f for f in file_names])
            self.print(f'{files}')

            self.counter = Counter(0)
            self.start_s = time.perf_counter()

        def finish_sample(self):
            elapsed_s = time.perf_counter() - self.start_s
            self.print(f'    {self.term.italic}Finshed: {self.term.normal}{elapsed_s:0.4f}s')
            del self.counter

        def update(self, t, sequence_t, n_passed):
            term = self.term
            if self.counter != 0:
                 self.print(term.move_up * 3, term.clear_eos, term.move_x(0), end='')
            self.print(f'       {term.bold}sequence:          {term.normal}{sequence_t:,}')
            self.print(f'       {term.bold}k-mers:            {term.normal}{t:,}')
            self.print(f'       {term.bold}sequences passed:  {term.normal}{n_passed:,}', end='')
            self.print(f' ({n_passed / sequence_t * 100.0:2f}%)')
            self.counter += 1

    def execute(self, args):
        for sample, name in iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names):
            self.status.start_sample(name, sample)
            for n_seqs, time, n_skipped in self.processor.chunked_process(*sample):
                self.status.update(time, n_seqs, self.processor.n_passed())
            self.status.finish_sample()

    def teardown(self):
        pass

