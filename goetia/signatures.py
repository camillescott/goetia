#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 15.10.2019

import os
import sys

import blessings
import curio
from pyfiglet import Figlet

from sourmash import SourmashSignature, save_signatures

from goetia import libgoetia, __version__
from goetia.cli import CommandRunner, get_output_interval_args
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.processors import AsyncSequenceProcessor
from goetia.messages import (Interval, DistanceCalc, SampleStarted, SampleFinished,
                             SampleSaturated, Error)
from goetia.saturation import SaturationPolicies
from goetia.utils import Namespace

SourmashSketch = libgoetia.signatures.SourmashSignature
UnikmerSketch  = libgoetia.signatures.UnikmerSignature


class SourmashRunner(CommandRunner):

    def __init__(self, parser):
        group = get_fastx_args(parser)
        group.add_argument('-i', dest='inputs', nargs='+', required=True)
        parser.add_argument('-K', default=31, type=int)
        parser.add_argument('-N', default=2500, type=int)
        parser.add_argument('--scaled', default=0, type=int)

        parser.add_argument('--save-sig',
                            help='Save the final, saturated signature to '
                                 'the given filename.')
        parser.add_argument('--save-stream',
                            action='store_true',
                            default=False,
                            help='Save the entire stream of signatures to '
                                 'the filename given in save-sig.')

        parser.add_argument('--distance-output', nargs='?')
        parser.add_argument('--cutoff', default=0.999, type=float)
        parser.add_argument('--saturation-policy', choices=list(SaturationPolicies.keys()),
                            default='median.distance')

        parser.add_argument('--tick-length', default=100000, type=int,
                            help='Approx. number of k-mers in a tick.')
        parser.add_argument('--window-size', default=5, type=int,
                            help='Number of ticks in a window.')

        parser.add_argument('--echo', default=None,
                            help='echo all events to the given file.')
        parser.add_argument('--term-graph', default=False, action='store_true',
                            help='draw a live distance graph the terminal.')

        parser.add_argument('--curio-monitor', default=False, action='store_true',
                            help='Run curio kernel monitor for async debugging.')

        super().__init__(parser)

    def postprocess_args(self, args):
        if args.term_graph:
            self.term_graph = Namespace()
            self.term_graph.term = blessings.Terminal()
            # minhash sigs only use jaccard distance
            args.distance_metric = 'jaccard'
        if args.save_stream and not args.save_sig:
            print('--save-stream requires --save-sig', file=sys.stderr)
            sys.exit(1)

    def make_signature(self, args):
        if args.scaled:
            return SourmashSketch.Signature.build(0, args.K, False, False, False, 42, args.scaled)
        else:
            return SourmashSketch.Signature.build(args.N, args.K, False, False, False, 42, 0)

    def setup(self, args):
        # build the underlying Processor specialized for sourmash signature
        self.signature = self.make_signature(args)
        processor = SourmashSketch.Processor.build(self.signature,
                                                   1000,
                                                   10000,
                                                   100000)
        # get the sample iter
        sample_iter = iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names)

        # build and save the async sequence processor
        self.processor = AsyncSequenceProcessor(processor, sample_iter, args.echo)
        
        # set up the saturation tracker
        def dfunc(sig_a, sig_b):
            return sig_a.similarity(sig_b)
        self.tracker = SaturationPolicies[args.saturation_policy](args.tick_length,
                                                                  args.window_size,
                                                                  args.cutoff,
                                                                  dfunc)
        self.sigs = []

        # set up a callback from Interval events on the sequence processor
        def on_interval(msg, events_q, args, sigs):

            kmer_time = self.processor.processor.n_inserted()
            sig = self.signature.to_sourmash()
            distance, delta, stat = self.tracker.push(sig, msg.t, kmer_time)

            if distance is not None:
                out_msg = DistanceCalc(sample_name=msg.sample_name,
                                       t=msg.t,
                                       delta=delta,
                                       distance=distance,
                                       stat=stat,
                                       stat_type=args.saturation_policy,
                                       file_names=msg.file_names)
                events_q.put(out_msg)

            if args.save_stream and not self.tracker.saturated:
                sigs.append(SourmashSignature(sig,
                                              name=f'{msg.sample_name}:{msg.t}',
                                              filename=format_filenames(msg.file_names)))

            if self.tracker.saturated:
                events_q.put(SampleSaturated(t=msg.t,
                                             sample_name=msg.sample_name,
                                             file_names=msg.file_names))
                self.processor.saturate()


        def on_saturated(msg, events_q, sigs):
            sigs.append(SourmashSignature(self.signature.to_sourmash(),
                                          name=f'{msg.sample_name}:{msg.t}:saturated',
                                          filename=format_filenames(msg.file_names)))

        def on_stop(msg, events_q, sigs):
            sigs.append(SourmashSignature(self.signature.to_sourmash(),
                                          name=f'{msg.sample_name}:{msg.t}:unsaturated',
                                          filename=format_filenames(msg.file_names)))
        
        # set up the listener and add the callback
        self.worker_listener = self.processor.add_listener('worker_q', 'sourmash.consumer')
        self.worker_listener.on_message(Interval, on_interval,
                                        self.processor.events_q, args, self.sigs)
        self.worker_listener.on_message(Error, on_stop,
                                        self.processor.events_q, self.sigs)

        self.events_listener = self.processor.add_listener('events_q', 'sourmash.producer')
        self.events_listener.on_message(SampleFinished, on_stop,
                                        self.processor.events_q, self.sigs)

        self.events_listener.on_message(SampleSaturated, on_saturated,
                                        self.processor.events_q, self.sigs)


        
        if args.term_graph:
            self.init_term_graph(args)
    
    def execute(self, args):
        if args.term_graph:
            with self.term_graph.term.hidden_cursor():
                curio.run(self.processor.start, with_monitor=args.curio_monitor)
        else:
            curio.run(self.processor.start, with_monitor=args.curio_monitor)

        if args.save_sig:
            with open(args.save_sig, 'w') as fp:
                save_signatures(self.sigs, fp=fp)

    def teardown(self):
        pass

    def init_term_graph(self, args):
        self.term_graph.frame = SignatureStreamFrame(self.term_graph.term, args)
        frame = self.term_graph.frame

        # set up frame drawing callbacks
        async def on_samplestart(msg):
            draw = await curio.spawn(frame.draw, messages=[f'Started processing on {msg.samplename}'], draw_dist_plot=False)
            await draw.join()
        self.term_graph.ss_listener = self.processor.add_listener('events_q', 'framedraw.samplestart')
        self.term_graph.ss_listener.on_message(SampleStarted, on_samplestart)
        
        async def on_distancecalc(msg):
            draw = await curio.spawn(frame.draw, msg.t, msg.distance, msg.stat)
            await draw.join()
        self.term_graph.dc_listener = self.processor.add_listener('events_q', 'framedraw.distancecalc')
        self.term_graph.dc_listener.on_message(DistanceCalc, on_distancecalc)


def format_filenames(file_names):
    return '+'.join([os.path.basename(name) for name in file_names])



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


class SignatureStreamFrame:

    def __init__(self, term, args, component_name='sourmash stream'):
        try:
            import plotille as pt
        except ImportError:
            print('plotille is required for terminal plotting: '
                  'install with `pip install plotille`',
                  file=sys.stderr)
            sys.exit(1)

        self.term = term
        self.args   = args

        self.distance_figure = pt.Figure()
        self.distance_figure.width = 60
        self.distance_figure.height = 20
        self.distance_figure.x_label = '# sequences'
        self.distance_figure.y_label = 'Δdistance'
        self.distance_figure.set_x_limits(min_=0)
        self.distance_figure.set_y_limits(min_=0)
        self.distance_figure.color_mode = 'byte'
        self.distance_figure.origin = False

        self.hist_figure = pt.Figure()
        self.hist_figure.width = 30
        self.hist_figure.height = 20
        self.hist_figure.x_label = 'Δdistance'
        self.hist_figure.y_label = 'counts'
        self.hist_figure.color_mode = 'byte'

        name_block = Figlet(font='doom').renderText('goetia').rstrip().split('\n')
        name_block[-1] = name_block[-1] + f' {__version__}'
        name_block.append(f'\n{term.bold}{component_name}{term.normal}')
        self.name_block = TextBlock('\n'.join(name_block))

        param_block = term.normal + term.underline + (' ' * 40) + term.normal + '\n\n'
        param_block += f'{term.bold}window size:       {term.normal}{args.window_size}\n'\
                       f'{term.bold}tick length:       {term.normal}{args.tick_length}\n'\
                       f'{term.bold}distance metric:   {term.normal}{args.distance_metric}\n'\
                       f'{term.bold}saturation policy: {term.normal}{args.saturation_policy}\n'\
                       f'{term.bold}cutoff:            {term.normal}{args.cutoff}'

        self.param_block = TextBlock(param_block)
        metric_text = term.normal + term.underline + (' ' * 40) + term.normal + '\n'
        metric_text += '{term.move_down}{term.bold}'
        metric_text += '# sequences:  '.ljust(20)
        metric_text += '{term.normal}{n_reads:,}\n'
        metric_text += '{term.bold}'
        metric_text += 'Δdistance:    '.ljust(20)
        metric_text += '{term.normal}{prev_d:.6f}\n'
        metric_text += '{term.bold}'
        metric_text += 'windowed(Δdistance): '.ljust(20)
        metric_text += '{term.normal}{stat:.6f}\n'
        metric_text += '{term.underline}' + (' ' * 40) + '{term.normal}'
        self.metric_text = metric_text

        #message_text = term.normal + term.underline + (' ' * 40) + term.normal + '\n'
        message_text = '{messages}'
        self.message_text = message_text

        self.distances = []
        self.distances_t = []

    def __del__(self):
        self._print(self.term.move_down * (self.distance_figure.height + 10))

    def _print(self, *args, **kwargs):
        self.term.stream.write(*args, **kwargs)

    def metric_block(self, n_reads, prev_d, stat):
        return TextBlock(self.metric_text.format(n_reads=n_reads, prev_d=prev_d, stat=stat, term=self.term))

    def message_block(self, messages):
        text = self.message_text.format(messages='\n'.join(messages))
        return TextBlock(text)

    def figure_block(self, distances, distances_t):
        self.distance_figure.clear()
        self.distance_figure.plot(distances_t, distances, lc=10, interp=None)
        text = self.distance_figure.show()
        return TextBlock(text)

    def hist_block(self, distances):
        self.hist_figure.clear()
        distances = np.log(distances)
        self.hist_figure.histogram(distances, lc=20)
        text = self.hist_figure.show()
        return TextBlock(text)

    async def draw(self, n_reads=0, prev_d=1.0, stat=1.0,
                   messages=None, draw_dist_plot=True,
                   draw_dist_hist=False):

        term = self.term
        self.distances.append(prev_d)
        self.distances_t.append(n_reads)

        figure_thr = await curio.spawn_thread(self.figure_block, self.distances, self.distances_t)
        metrics_thr = await curio.spawn_thread(self.metric_block, n_reads, prev_d, stat)
        
        figure = await figure_thr.join()
        metrics = await metrics_thr.join()

        term.stream.write(term.clear_eos)

        if draw_dist_plot:
            with term.location():
                figure.draw(term, max(self.name_block.width, self.param_block.width) + 1, 1)

        with term.location():
            self.name_block.draw(term, 0, shift_y=1)

        with term.location():
            self.param_block.draw(term, 0, shift_y=self.name_block.height+4)

        with term.location():
            metrics.draw(term,
                         0,
                         self.name_block.height + self.param_block.height + 4)

        if messages is not None:
            with term.location():
                ypos = _metric_block.height + self.name_block.height + self.param_block.height + 9
                self.message_block(messages).draw(term, 0, ypos)

        if draw_dist_hist and distances is not None:
            with term.location():
                xpos = max(self.name_block.width, self.param_block.width) + 2
                if draw_dist_plot:
                    xpos += self.distance_figure.width + 5
                self.hist_block(distances).draw(term, xpos, 1)
