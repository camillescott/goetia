#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 15.10.2019

import curio

from boink import libboink
from boink.cli import CommandRunner, get_output_interval_args
from boink.parsing import get_pairing_args, iter_fastx_inputs
from boink.processors import AsyncSequenceProcessor
from boink.messages import Interval, DistanceCalc

SourmashSignature = libboink.signatures.SourmashSignature
UnikmerSignature  = libboink.signatures.UnikmerSignature


class SourmashRunner(CommandRunner):

    def __init__(self, parser):
        get_output_interval_args(parser)
        group = get_pairing_args(parser)
        group.add_argument('-i', dest='inputs', nargs='+', required=True)
        parser.add_argument('-K', default=31)
        parser.add_argument('-N', default=10000)
        super().__init__(parser)

    def postprocess_args(self, args):
        pass

    def setup(self, args):
        # build the sourmash signature
        self.signature = SourmashSignature.Signature.build(args.N, args.K, False, 42, 0)

        # build the underlying Processor specialized for sourmash signature
        processor = SourmashSignature.Processor.build(self.signature,
                                                      args.fine_interval,
                                                      args.medium_interval,
                                                      args.coarse_interval)

        # get the sample iter
        sample_iter = iter_fastx_inputs(args.inputs, args.pairing_mode)

        # build and save the async sequence processor
        self.processor = AsyncSequenceProcessor(processor, sample_iter)

        # set up a callback from Interval events on the sequence processor
        self.distances = []
        self.times = []
        self.prev_mh = {}

        def on_interval(msg, events_q):
            curr_mh = self.signature.to_sourmash()
            if self.prev_mh:
                self.distances.append(self.prev_mh.similarity(curr_mh))
                self.times.append(msg.t)
                delta = self.times[-1] - self.times[-2] if len(self.times) > 1 else self.times[-1]
                out_msg = DistanceCalc(sample_name=msg.sample_name,
                                       t=msg.t,
                                       delta=delta,
                                       distance=self.distances[-1])
                # note events_q is a UniversalQueue so doesn't need to be awaited
                events_q.put(out_msg)
            self.prev_mh = curr_mh
        
        # set up the listener and add the callback
        self.interval_listener = self.processor.add_listener('worker_q', 'distances')
        self.interval_listener.on_message(Interval, on_interval,
                                          self.processor.events_q)
    
    def execute(self, args):
        curio.run(self.processor.run, with_monitor=True)
    
    def teardown(self):
        pass
