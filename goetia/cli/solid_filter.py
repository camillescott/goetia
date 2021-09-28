#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : solid_filter.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 12.03.2020

import sys

from goetia.filters import SolidFilter
from goetia.dbg import get_graph_args, process_graph_args
from goetia.cli.args import get_output_interval_args
from goetia.cli.runner import CommandRunner
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.storage import get_storage_args, process_storage_args


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

    def execute(self, args):
        for sample, name in iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names):
            for n_seqs, time, n_skipped in self.processor.chunked_process(*sample):
                print(f'{sample}, {name}: {n_seqs} reads, {self.processor.n_passed()} passed filter.',
                      file=sys.stderr)

    def teardown(self):
        pass

