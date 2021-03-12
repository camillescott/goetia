#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 28.07.2020

import os
import sys

import blessings
import curio
import numpy as np

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
from goetia.utils import Namespace


class SignatureRunner(CommandRunner):

    def __init__(self, parser):
        get_output_interval_args(parser)
        group = get_fastx_args(parser)
        group.add_argument('-i', dest='inputs', nargs='+', required=True)
        #parser.add_argument('-K', default=31, type=int)
        #parser.add_argument('-N', default=2500, type=int)
        #parser.add_argument('--scaled', default=0, type=int)

        parser.add_argument('--save-sig',
                            help='Save the final, saturated signature to '
                                 'the given filename.')
        parser.add_argument('--save-stream',
                            action='store_true',
                            default=False,
                            help='Save the entire stream of signatures to '
                                 'the filename given in save-sig.')
        parser.add_argument('--save-stream-tick',
                           type=int,
                           default=10,
                           help='Save a copy of the signature every N intervals.')

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

    def _distance_func(self, sigs):
        """ Compute the distance between two signatures in sigs.

        Args:
            sigs (list): The signatures.
        
        Returns:
            float: Distance between the signatures.
        """
        raise NotImplementedError()

    @staticmethod
    def _make_signature(args):
        raise NotImplementedError()

    @staticmethod
    def _make_processor(signature, args):
        raise NotImplementedError()

    @staticmethod
    def _convert_signature(sig, msg):
        return sig
    
    @staticmethod
    def _save_signatures(sigs, args):
        raise NotImplementedError()

    @staticmethod
    def _on_interval(msg, events_q, args, runner):
        sig = runner._convert_signature(runner.signature, msg)
        distance, _ = runner.sigs.push((sig, msg.t))

        if not np.isnan(distance):
            cutoff_reached, stat, _ = runner.cutoff.push((distance, msg.t))

            if not np.isnan(stat):
                out_msg = DistanceCalc(sample_name=msg.sample_name,
                                       t=msg.t,
                                       sequence=msg.sequence,
                                       distance=distance,
                                       stat=stat,
                                       stat_type=runner.cutoff.name,
                                       file_names=msg.file_names)
                events_q.put(out_msg)

            if cutoff_reached and args.saturate:
                events_q.put(SampleSaturated(t=msg.t,
                                             sequence=msg.sequence,
                                             sample_name=msg.sample_name,
                                             file_names=msg.file_names))
                runner.processor.saturate()
    
    def setup(self, args):
        # Create the signature and the libgoetia sequence processor: implemented by subclass
        self.signature = self._make_signature(args)
        signature_processor = self._make_processor(self.signature, args)

        # Pass the libgoetia processor over the async handler
        self.processor = AsyncSequenceProcessor(signature_processor,
                                                iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names),
                                                echo=args.echo)
        
        # Signatures: held in a RollingPairwise, stores the signatures themselves (if desired)
        # and calls the distance function
        self.sigs = RollingPairwise(self._distance_func, history = int(not args.save_stream))
        
        # Smooths the distances over a window and checks the smooth values against
        # a cutoff function to trigger saturation
        self.cutoff = SlidingCutoff(args.window_size,
                                    args.smoothing_function_func,
                                    args.cutoff_function_func)

        #
        # Set up the event handlers. We set up two listeners:
        #     - The worker listener: registered on the worker queue of the async processor;
        #       this is the producer loaded by the libgoetia processor loop.
        #     - The events listener: this is loaded by other handlers. It has everything fom the worker
        #       queue and any events triggered by its subscribers.
        #

        self.worker_listener = self.processor.add_listener('worker_q', 'signature.consumer')
        self.worker_listener.on_message(Interval, self._on_interval,
                                        self.processor.events_q, args, self)

        if args.term_graph:
            self.init_term_graph_io(args)
        else:
            self.init_std_io(args)
    
    def execute(self, args):
        if args.term_graph:
            with self.term_graph.term.hidden_cursor():
                curio.run(self.processor.start, with_monitor=args.curio_monitor)
        else:
            curio.run(self.processor.start, with_monitor=args.curio_monitor)

        if args.save_sig:
            if args.save_stream:
                self._save_signatures(self.sigs.values(), args)
            else:
                self._save_signatures([self.sigs.tail()], args)

    def teardown(self):
        pass

    def init_term_graph_io(self, args):
        self.term_graph.frame = SignatureStreamFrame(self.term_graph.term, args)
        frame = self.term_graph.frame

        self.term_graph.listener = self.processor.add_listener('events_q', 'signature.term_graph')

        # set up frame drawing callbacks
        async def on_samplestart(msg):
            draw = await curio.spawn(frame.draw, 0, 0, 1.0, 1.0, [f'Started processing on {msg.sample_name}'], False)
            await draw.join()

        async def on_distancecalc(msg):
            draw = await curio.spawn(frame.draw, msg.t, msg.sequence, msg.distance, msg.stat)
            await draw.join()

        self.term_graph.listener.on_message(SampleStarted, on_samplestart)
        self.term_graph.listener.on_message(DistanceCalc, on_distancecalc)

    @staticmethod
    def _on_saturated(msg):
        print(f'Signature saturated: {msg.sample_name}->{msg.file_names} at time={msg.t}, sequence={msg.sequence}', file=sys.stderr)

    @staticmethod
    def _on_error(msg):
        print(f'ERROR: {msg.sample_name}->{msg.file_names} at time={msg.t}, sequence={msg.sequence} '\
              f'\n-- BEGIN EXCEPTION --\n{msg.error}\n-- END EXCEPTION --', file=sys.stderr)

    @staticmethod
    def _on_sample_finished(msg):
        print(f'Sample finished: {msg.sample_name}->{msg.file_names} at time={msg.t}, sequence={msg.sequence}', file=sys.stderr)
    
    @staticmethod
    def _on_status_update(msg):
        print(f'{msg.sample_name}: distance={msg.distance} at sequence={msg.sequence} and time={msg.t}...', file=sys.stderr)
    
    def init_std_io(self, args):
        self.events_listener = self.processor.add_listener('events_q', 'signature.std_io')

        self.events_listener.on_message(Error, self._on_error)
        self.events_listener.on_message(SampleFinished, self._on_sample_finished)
        self.events_listener.on_message(SampleSaturated, self._on_saturated)
        self.events_listener.on_message(DistanceCalc, self._on_status_update)
