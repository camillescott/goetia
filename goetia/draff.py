#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : draff.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import functools
import os
import sys

from pyfiglet import Figlet
import scipy
import numpy as np
from scipy.spatial import distance as dmetrics

from cppyy.gbl import std

from goetia import libgoetia
from goetia.hashing import Canonical, StrandAware
from goetia.signatures import DraffSignature
from goetia.storage import get_storage_args, process_storage_args
from goetia.cli.signature_runner import SignatureRunner


class DraffStreamFrame:

    def __init__(self, term, args):
        self.args   = args

        self.distance_figure = pt.Figure()
        self.distance_figure.width = 60
        self.distance_figure.height = 20
        self.distance_figure.x_label = '# sequences'
        self.distance_figure.y_label = 'Δdistance'
        self.distance_figure.set_x_limits(min_=0)
        self.distance_figure.set_y_limits(min_=0)
        self.distance_figure.color_mode = 'byte'

        self.hist_figure = pt.Figure()
        self.hist_figure.width = 30
        self.hist_figure.height = 20
        self.hist_figure.x_label = 'Δdistance'
        self.hist_figure.y_label = 'counts'
        self.hist_figure.color_mode = 'byte'

        self.term = term

        name_block = Figlet(font='doom').renderText('draff').rstrip().split('\n')
        name_block[-1] = name_block[-1] + ' stream, v0.1'
        self.name_block = TextBlock('\n'.join(name_block))

        param_block = term.normal + term.underline + (' ' * 40) + term.normal + '\n\n'
        param_block += '{term.bold}distance window size:  {term.normal}{distance_window}\n'\
                       '{term.bold}distance metric:       {term.normal}{distance_metric}\n'\
                       '{term.bold}distance stdev cutoff: {term.normal}{stdev_cutoff}'.format(term = self.term,
                                                                                              distance_window = args.distance_window,
                                                                                              stdev_cutoff    = args.stdev_cutoff)
class DraffRunner(SignatureRunner):

    def __init__(self, parser):
        get_storage_args(parser)
        parser.add_argument('-W', type=int, default=31)
        parser.add_argument('-K', type=int, default=9)
        parser.add_argument('--prefix', nargs='*')
        parser.add_argument('--merge')
        parser.add_argument('--distance-metric', default='cosine',
                            choices=['cosine', 'euclidean', 'canberra', 'braycurtis',
                                     'chebyshev', 'correlation', 'minkowski',
                                     'cityblock', 'jensenshannon', 'mahalanobis',
                                     'seuclidean', 'sqeuclidean'])

        super().__init__(parser)
    
    def postprocess_args(self, args):
        process_storage_args(args)
        args.signature_t = libgoetia.signatures.UnikmerSignature[args.storage, StrandAware]
        super().postprocess_args(args)
        self._distance_metric = getattr(dmetrics, args.distance_metric)
    
    def _distance_func(self, sigs):
        sig_a, sig_b = sigs
        return self._distance_metric(sig_a.sketch, sig_b.sketch)

    @staticmethod
    def _make_signature(args):
        return args.signature_t.Signature.build(args.W, args.K, storage_args=args.storage_args)
    
    @staticmethod
    def _make_processor(signature, args):
        return args.signature_t.Processor.build(signature,
                                                args.interval)
    
    @staticmethod
    def _convert_signature(sig, msg):
        return DraffSignature(sig, name=msg.sample_name)

    @staticmethod
    def _save_signatures(sigs, args):

        with open(prefix + '.draffsig.json', 'w') as fp:
            sig_gen.save(fp, prefix)

    def save_distances(self, prefix, distances, distances_t):
        with open(prefix + '.distances.csv', 'w') as fp:
            print('t', ',', 'distance', file=fp)
            for t, distance in zip(distances_t, distances):
                print(t, ',', distance, file=fp)


def run_draff():
    term   = blessings.Terminal()
    with term.hidden_cursor():
        runner = DraffRunner()
        stream = DraffStream(runner, term)

        runner.run()

    return 0
