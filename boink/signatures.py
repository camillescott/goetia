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
from boink.processors import AsyncSequenceProcessor, MessageTypes

SourmashSignature = libboink.signatures.SourmashSignature
UnikmerSignature  = libboink.signatures.UnikmerSignature



class NewDistanceCalc:
    @staticmethod
    def generate(name, t, delta, distance):
        assert t >= 0
        return 'NewDistanceCalc', {'t': t, 'delta': delta, 'samople.name': name, 'distance': distance}


class AsyncSourmashProcessor(AsyncSequenceProcessor):

    def __init__(self, args, N=10000, K=31):
        self.signature = SourmashSignature.Signature.build(N, K, False, 42, 0)
        processor = SourmashSignature.Processor.build(self.signature,
                                                      args.fine_interval,
                                                      args.medium_interval,
                                                      args.coarse_interval)
        sample_iter = iter_fastx_inputs(args.inputs, args.pairing_mode)
        super().__init__(processor, sample_iter)

    async def track_distances(self):
        distances = []
        times = []
        prev_mh = None

        try:
            msg_q = curio.Queue()
            self.subscribe('worker_q', msg_q, 'track_distances')
            while True:
                msg = await msg_q.get()

                if msg is None:
                    await msg_q.task_done()
                    break

                if msg['MSG_TYPE'] == 'Interval':
                    curr_mh = self.signature.to_sourmash()
                    if prev_mh is not None:
                        distances.append(prev_mh.similarity(curr_mh))
                        times.append(msg['DATA']['t'])
                        delta = times[-1] - times[-2] if len(times) > 1 else times[-1]
                        await self.events_q.put(MessageTypes.generate(NewDistanceCalc,
                                                                      msg['DATA']['sample.name'],
                                                                      msg['DATA']['t'],
                                                                      delta,
                                                                      distances[-1]))
                    prev_mh = curr_mh

                await msg_q.task_done()

        except curio.CancelledError:
            raise
        finally:
            self.unsubscribe('events_q', msg_q)


    async def run(self):
        await super().run(extra_tasks=[self.track_distances])
        #await super().run()


class SourmashRunner(CommandRunner):

    def __init__(self, parser):
        get_output_interval_args(parser)
        group = get_pairing_args(parser)
        group.add_argument('-i', dest='inputs', nargs='+', required=True)
        super().__init__(parser)

    def postprocess_args(self, args):
        pass

    def setup(self, args):
        self.async_proc = AsyncSourmashProcessor(args)
    
    def execute(self, args):
        curio.run(self.async_proc.run, with_monitor=True)
    
    def teardown(self):
        pass
