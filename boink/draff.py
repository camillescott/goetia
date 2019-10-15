#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : draff.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import argparse
import functools
import os
import sys

import blessings
from pyfiglet import Figlet
import plotille as pt
import scipy
import numpy as np
from scipy.spatial import distance as dmetrics

from cppyy.gbl import std

from boink.cli import get_output_interval_args 
from boink.parsing import get_pairing_args, iter_fastx_inputs
from boink.data import load_unikmer_map
from boink import libboink
from boink.storage import get_storage_args, process_storage_args
from boink.utils import find_common_basename, remove_fx_suffix, grouper


class TextBlock:

    def __init__(self, text):
        self.text = text.split('\n')
        self.height = len(self.text)
        self.width = max((len(line) for line in self.text))

    def __iter__(self):
        yield from self.text

    def draw(self, term, x, shift_y=None):
        out = term.stream
        
        if shift_y:
            buf = term.move_down * shift_y
        else:
            buf = ''

        for line in self.text:
            buf += term.move_x(x)
            buf += line
            buf += term.move_down
        
        out.write(buf)


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
                                                                                              distance_metric = args.distance_metric,
                                                                                              stdev_cutoff    = args.stdev_cutoff)
        self.param_block = TextBlock(param_block)
        metric_text = term.normal + term.underline + (' ' * 40) + term.normal + '\n'
        metric_text += '{term.move_down}{term.bold}'
        metric_text += '# sequences:  '.ljust(20)
        metric_text += '{term.normal}{n_reads:,}\n'
        metric_text += '{term.bold}'
        metric_text += 'Δdistance:    '.ljust(20)
        metric_text += '{term.normal}{prev_d:.20f}\n'
        metric_text += '{term.bold}'
        metric_text += 'σ(Δdistance): '.ljust(20)
        metric_text += '{term.normal}{stdev:.20f}\n'
        metric_text += '{term.underline}' + (' ' * 40) + '{term.normal}'
        self.metric_text = metric_text

        #message_text = term.normal + term.underline + (' ' * 40) + term.normal + '\n'
        message_text = '{messages}'
        self.message_text = message_text


    def __del__(self):
        self._print(self.term.move_down * (self.distance_figure.height + 10))

    def _print(self, *args, **kwargs):
        self.term.stream.write(*args, **kwargs)

    def metric_block(self, n_reads, prev_d, stdev):
        return TextBlock(self.metric_text.format(n_reads=n_reads, prev_d=prev_d, stdev=stdev, term=self.term))

    def message_block(self, messages):
        text = self.message_text.format(messages='\n'.join(messages))
        return TextBlock(text)

    def figure_block(self, distances, distances_t):
        self.distance_figure.clear()
        self.distance_figure.plot(distances_t, distances, lc=10)
        text = self.distance_figure.show()
        return TextBlock(text)

    def hist_block(self, distances):
        self.hist_figure.clear()
        distances = np.log(distances)
        self.hist_figure.histogram(distances, lc=20)
        text = self.hist_figure.show()
        return TextBlock(text)

    def draw(self, n_reads=0, prev_d=1.0, stdev=1.0, distances=None, distances_t=None,
                   messages=None, draw_dist_plot=True, draw_dist_hist=False):

        term = self.term
        term.stream.write(term.clear_eos)

        with term.location():
            self.name_block.draw(term, 0, shift_y=1)

        with term.location():
            self.param_block.draw(term, 0, shift_y=self.name_block.height+4)

        _metric_block = self.metric_block(n_reads, prev_d, stdev)
        with term.location():
            _metric_block.draw(term,
                               0,
                               self.name_block.height + self.param_block.height + 4)

        if messages is not None:
            with term.location():
                ypos = _metric_block.height + self.name_block.height + self.param_block.height + 9
                self.message_block(messages).draw(term, 0, ypos)


        if draw_dist_plot and None not in (distances, distances_t):
            with term.location():
                self.figure_block(distances, distances_t).draw(term, max(self.name_block.width, self.param_block.width) + 1, 1)

        if draw_dist_hist and distances is not None:
            with term.location():
                xpos = max(self.name_block.width, self.param_block.width) + 2
                if draw_dist_plot:
                    xpos += self.distance_figure.width + 5
                self.hist_block(distances).draw(term, xpos, 1)


class DraffRunner:

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('--version', action='version', version='draff-0.1')
        self.parser.set_defaults(func = lambda _: self.parser.print_help())
        self.subparsers = self.parser.add_subparsers()

    def add_subparser(self, name, func):
        parser = self.subparsers.add_parser(name)
        parser.add_argument('-W', type=int, default=31)
        parser.add_argument('-K', type=int, default=9)
        parser.add_argument('--prefix', nargs='*')
        get_storage_args(parser)
        parser.add_argument('--merge')
        parser.add_argument('-i', '--inputs', nargs='+', required=True)

        get_pairing_args(parser)
        get_output_interval_args(parser)
        parser.set_defaults(func=func)

        return parser

    def run(self):
        args = self.parser.parse_args()
        process_storage_args(args)
        self.signature_t = libboink.signatures.UnikmerSignature[args.storage]
        return args.func(args)

    def iter_inputs(self, args):
        if args.pairing_mode == 'split':
            _samples  = list(grouper(2, args.inputs))
            _prefixes = [find_common_basename(remove_fx_suffix(l),
                                              remove_fx_suffix(r)) for l, r in _samples]
        else:
            _samples  = [(s, ) for s in args.inputs]
            _prefixes = [os.path.basename(remove_fx_suffix(s)) for s, in _samples]

        if args.merge:
            _prefixes = [args.merge] * len(_samples)

        if args.prefix:
            if args.merge:
                pass
            else:
                if len(args.prefix) != len(_samples):
                    raise RuntimeError('Number of names must match number of samples.')
                _prefixes = args.prefix

        yield from zip(_samples, _prefixes)


class DraffStream:

    def __init__(self, draff_runner, term):
        self.runner = draff_runner
        self.parser = draff_runner.add_subparser('stream', self.run)
        self.parser.add_argument('--distance-metric', default='cosine')
        self.parser.add_argument('--distance-output', nargs='?')
        self.parser.add_argument('--stdev-cutoff', default=0.00001, type=float)
        self.parser.add_argument('--distance-window', default=5, type=int)

        self.term   = term
        self._print = functools.partial(print, end='')

    def run(self, args):
        frame = DraffStreamFrame(self.term, args)
        term = frame.term
        frame.draw()

        dfunc = getattr(dmetrics, args.distance_metric)

        frame.draw(messages=[term.read + 'Initializing UKHS({},{})...'.format(args.W, args.K)],
                   draw_dist_plot=False)
        ukhs = load_unikmer_map(args.W, args.K)
        signature_t = self.runner.signature_t
        
        if args.merge:
            frame.draw(messages=[term.red + 'Initializing signature for {}...'.format(args.merge)], 
                       draw_dist_plot=False)
            sig_gen = signature_t.Signature.build(args.W, args.K, ukhs, *args.storage_args)
            processor = signature_t.Processor.build(sig_gen,
                                                    args.fine_interval,
                                                    args.medium_interval,
                                                    args.coarse_interval)

            for sample, prefix in self.runner.iter_inputs(args):
                frame.draw(messages=['Started processing on {}'.format(sample)], draw_dist_plot=False)

                n_reads, distances, distances_t, stdevs, saturated = self.process_sample(sig_gen, processor, sample, frame,
                                                                                        dfunc, args.distance_window, args.stdev_cutoff)
                if saturated:
                    msgs = [term.green
                            + '{name} reached saturation at read {n:,} in sample {sample}.'.format(name=prefix, 
                                                                                                   n=n_reads,
                                                                                                   sample=sample)
                            + term.normal]
                else:
                    msgs = [term.red   + '{name} didn\'t saturate at {sample}; {n:,} reads in sample.'.format(name=name, 
                                                                                                              sample=sample,
                                                                                                              n=n_reads) + term.normal]
                frame.draw(n_reads, distances[-1], stdevs[-1], distances, distances_t, messages=msgs)

                if saturated:
                    break

            self.save_distances(args.merge, distances, distances_t)
            args.save_signature(args.merge, sig_gen)

        else:
            for sample, prefix in self.runner.iter_inputs(args):
                frame.draw(messages=[term.red + 'Initializing signature for {}...'.format(prefix)], 
                           draw_dist_plot=False)
                sig_gen = signature_t.Signature.build(args.W, args.K, ukhs.__smartptr__(), *args.storage_args)
                processor = signature_t.Processor.build(sig_gen,
                                                        args.fine_interval,
                                                        args.medium_interval,
                                                        args.coarse_interval)

                frame.draw(messages=['Started processing on {}.'.format(sample)])
                n_reads, distances, distances_t, stdevs,  saturated = self.process_sample(sig_gen, processor, sample, frame,
                                                                                          dfunc, args.distance_window, args.stdev_cutoff)

                if saturated:
                    msgs = [term.green + '{name} reached saturation at {n:,} reads.'.format(name=prefix, n=n_reads) + term.normal]
                else:
                    msgs = [term.red   + '{name} didn\'t saturate; {n:,} reads in sample.'.format(name=prefix, n=n_reads) + term.normal]

                frame.draw(n_reads, distances[-1], stdevs[-1], distances, distances_t, messages=msgs)

                self.save_distances(prefix, distances, distances_t)
                self.save_signature(prefix, sig_gen)

        self._print(term.move_down)

    def process_sample(self, sig_gen, processor, sample, frame, dfunc, distance_window, stdev_cutoff):
        saturated   = False
        signatures  = []
        distances   = []
        distances_t = []
        stdevs      = []

        for n_reads, state in processor.chunked_process(*sample):
            if state.fine:
                signatures.append(sig_gen.signature)
                frame.draw(n_reads, messages=['Accumulating distances over initial window.'])

                if len(signatures) > 1:
                    distances.append(dfunc(signatures[-2], signatures[-1]))
                    distances_t.append(n_reads)
                
                if len(distances) >= distance_window:
                    stdev = np.std(distances[-distance_window:])
                    stdevs.append(stdev)
                    if len(stdevs) >= distance_window and \
                            all((d < stdev_cutoff for d in stdevs[-distance_window:])):
                        
                       return n_reads, distances, distances_t, stdevs, True

                    frame.draw(n_reads, distances[-1], stdevs[-1], distances, distances_t)

        return n_reads, distances, distances_t, stdevs, False

    def save_distances(self, prefix, distances, distances_t):
        with open(prefix + '.distances.csv', 'w') as fp:
            print('t', ',', 'distance', file=fp)
            for t, distance in zip(distances_t, distances):
                print(t, ',', distance, file=fp)

    def save_signature(self, prefix, sig_gen):
        with open(prefix + '.draffsig.json', 'w') as fp:
            sig_gen.save(fp, prefix)


def run_draff():
    term   = blessings.Terminal()
    with term.hidden_cursor():
        runner = DraffRunner()
        stream = DraffStream(runner, term)

        runner.run()

    return 0
