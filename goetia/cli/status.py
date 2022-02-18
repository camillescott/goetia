#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : status.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 17.02.2022

import sys

import blessings

from goetia.utils import Counter


class StatusOutput:

    def __init__(self, names, term=None, file=sys.stderr):
        self.names = names
        self.file = file
        self.term = term or blessings.Terminal()

    def print(self, *args, **kwargs):
        print(*args, **kwargs, file=self.file)

    def message(self, *args, **kwargs):
        self.counter = Counter(0)
        self.print(*args, **kwargs)

    def start_sample(self, sample_name, file_names):
        self.print(f'{self.term.italic}Begin sample: {self.term.normal}{sample_name}')
        files = '\n'.join(['    ' + f for f in file_names])
        self.print(f'{files}')

        self.counter = Counter(0)
        self.interval_start_t = 0
        self.interval_start_seq = 0

    def finish_sample(self, seconds_elapsed):
        self.print(f'    {self.term.italic}Finished: {self.term.normal}{seconds_elapsed:0.4f}s')
        del self.counter

    def update(self, t, sequence_t, seconds_elapsed, *args):
        if not hasattr(self, 'counter'):
            return
        if len(args) != len(self.names):
            raise ValueError(f'Got {len(args)} values, expected {len(self.names)}.')

        interval_t = t - self.interval_start_t
        interval_kmers_per_s = int(interval_t / seconds_elapsed)
        interval_seqs = sequence_t - self.interval_start_seq
        interval_seqs_per_s = int(interval_seqs / seconds_elapsed)

        names = ['sequence: ', 'k-mers: ', 'kmers/s: ', 'seqs/s: ']
        for name in self.names:
            names.append(f'{name}: ')

        values = [f'{sequence_t:,}', f'{t:,}', f'{interval_kmers_per_s:,}', f'{interval_seqs_per_s:,}']
        for value in args:
            if isinstance(value, str):
                values.append(value)
            elif isinstance(value, float):
                values.append(f'{value:0.7f}')
            elif isinstance(value, int):
                values.append(f'{value:,}')
            else:
                values.append(str(value))

        term = self.term
        if self.counter != 0:
             self.print(term.move_up * len(values), term.clear_eos, term.move_x(0), end='')
        width = max(map(len, names))

        for name, value in zip(names, values):
            self.print(f'        {term.bold}{name:<{width}}{term.normal}{value}')

        self.counter += 1
        self.interval_start_t = t
        self.interval_start_seq = sequence_t
