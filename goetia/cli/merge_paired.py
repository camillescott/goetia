#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : merge_paired.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 29.10.2021

import sys
import time

from goetia import libgoetia
from goetia.alphabets import DNA_SIMPLE, DNAN_SIMPLE
from goetia.cli.runner import CommandRunner
from goetia.cli.args import get_output_interval_args
from goetia.parsing import FastxParser, get_fastx_args, iter_fastx_inputs
from goetia.utils import Counter

import blessings


desc = '''
{term.italic}utils module: {term.normal}merge paired

    Merge split reads into a single file.
'''


class PairedMergeRunner(CommandRunner):

    def __init__(self, parser):
        get_output_interval_args(parser)
        group = get_fastx_args(parser, default_mode='split')
        group.add_argument('-o', dest='output_filename', default='/dev/stdout')
        group.add_argument('-i', '--inputs', dest='inputs', nargs='+', required=True)
        group.add_argument('--quiet', action='store_true', default=False)

        super().__init__(parser, description=desc)

    def postprocess_args(self, args):
        pass

    def setup(self, args):
        self.processor = libgoetia.MergeProcessor[FastxParser[DNAN_SIMPLE]].build(args.output_filename,
                                                        args.interval)
        self.status = PairedMergeRunner.StatusOutput(term=self.term, quiet=args.quiet)

    def execute(self, args):
        for sample, name in iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names):
            self.status.start_sample(name, sample)
            for n_seqs, time, n_skipped in self.processor.chunked_process(*sample):
                self.status.update(n_seqs)
            self.status.finish_sample()

    def teardown(self):
        pass

    class StatusOutput:

        def __init__(self, term=None, file=sys.stderr, quiet=False):
            self.file = file
            self.term = term or blessings.Terminal()
            self.quiet = quiet

        def print(self, *args, **kwargs):
            if not self.quiet:
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

        def update(self, sequence_t):
            term = self.term
            if self.counter != 0:
                 self.print(term.move_up, term.clear_eos, term.move_x(0), end='')
            self.print(f'       {term.bold}sequence:          {term.normal}{sequence_t:,}')
            self.counter += 1
