#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : sourmash_stream.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 04.06.2020

import os
import sys

import blessings
import curio
import numpy as np
from pyfiglet import Figlet

from sourmash import SourmashSignature, save_signatures

from goetia import libgoetia, __version__
from goetia.cli.runner import CommandRunner
from goetia.cli.args import get_output_interval_args
from goetia.cli.cli import format_filenames, TextBlock
from goetia.cli.signature_frame import SignatureStreamFrame
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.processors import AsyncSequenceProcessor
from goetia.messages import (Interval, DistanceCalc, SampleStarted, SampleFinished,
                             SampleSaturated, Error)
from goetia.saturation import (cutoff_functions, smoothing_functions,
                               SlidingCutoff, RollingPairwise)
from goetia.signatures import SourmashSketch
from goetia.utils import Namespace


class SourmashRunner(CommandRunner):

    def __init__(self, parser):
        get_output_interval_args(parser)
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

        parser.add_argument('--smoothing-function',
                            choices=list(smoothing_functions.keys()),
                            default='mean')
        parser.add_argument('--cutoff-function',
                            choices=list(cutoff_functions.keys()),
                            default='all')
        parser.add_argument('--cutoff', default=0.999, type=float)
        parser.add_argument('--saturate', default=False, action='store_true')
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
        args.smoothing_function_func = smoothing_functions[args.smoothing_function]
        args.cutoff_function_func = cutoff_functions[args.cutoff_function](args.cutoff)

    def make_signature(self, args):
        if args.scaled:
            return SourmashSketch.Signature.build(0, args.K, False, False, False, 42, args.scaled)
        else:
            return SourmashSketch.Signature.build(args.N, args.K, False, False, False, 42, 0)

    def setup(self, args):
        # build the underlying Processor specialized for sourmash signature
        print('RUN SETUP')
        self.signature = self.make_signature(args)
        processor = SourmashSketch.Processor.build(self.signature,
                                                   args.interval)
        
        # get the sample iter
        sample_iter = iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names)

        # build and save the async sequence processor
        self.processor = AsyncSequenceProcessor(processor, sample_iter, args.echo)
        
        # set up the saturation tracker
        def dfunc(sigs):
            sig_a, sig_b = sigs
            sim = sig_a.similarity(sig_b)
            print('sim:', sim)
            return sim
        
        self.sigs = RollingPairwise(dfunc, history = int(not args.save_stream))
        
        self.cutoff = SlidingCutoff(args.window_size,
                                    args.smoothing_function_func,
                                    args.cutoff_function_func)

        # set up a callback from Interval events on the sequence processor
        def on_interval(msg, events_q, args, sigs):

            sig = SourmashSignature(self.signature.to_sourmash(),
                                    name=f'{msg.sample_name}:{msg.t}',
                                    filename=format_filenames(msg.file_names))
            distance, _ = sigs.push((sig, msg.t))

            if not np.isnan(distance):
                cutoff_reached, stat, _ = self.cutoff.push((distance, msg.t))

                if not np.isnan(stat):
                    out_msg = DistanceCalc(sample_name=msg.sample_name,
                                           t=msg.t,
                                           distance=distance,
                                           stat=stat,
                                           stat_type=self.cutoff.name,
                                           file_names=msg.file_names)
                    events_q.put(out_msg)

                if cutoff_reached and args.saturate:
                    events_q.put(SampleSaturated(t=msg.t,
                                                sample_name=msg.sample_name,
                                                file_names=msg.file_names))
                    self.processor.saturate()


        def on_saturated(msg, events_q, sigs):
            self.sigs.push((SourmashSignature(self.signature.to_sourmash(),
                                          name=f'{msg.sample_name}:{msg.t}:saturated',
                                          filename=format_filenames(msg.file_names)),
                            msg.t))

        def on_stop(msg, events_q, sigs):
            self.sigs.push((SourmashSignature(self.signature.to_sourmash(),
                                          name=f'{msg.sample_name}:{msg.t}:stopped',
                                          filename=format_filenames(msg.file_names)),
                            msg.t))
        
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
                save_signatures(self.sigs.values(), fp=fp)

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
