#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : streamhasher.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 17.02.2022

import curio

from goetia import libgoetia
from goetia.cli.runner import CommandRunner
from goetia.cli.args import get_output_interval_args
from goetia.cli.status import StatusOutput
from goetia.dbg import get_graph_args, process_graph_args
from goetia.messages import (Interval,
                             SampleStarted,
                             SampleFinished,
                             Error)
from goetia.parsing import get_fastx_args, iter_fastx_inputs
from goetia.processors import (AsyncSequenceProcessor,
                               AsyncJSONStreamWriter)


desc = '''
{term.italic}utils module: {term.normal}stream hasher
'''


class StreamHasherRunner(CommandRunner):
    
    def __init__(self, parser):
        get_output_interval_args(parser)
        get_graph_args(parser)
        group = get_fastx_args(parser)
        group.add_argument('-i', dest='inputs', nargs='+', required=True)

        parser.add_argument('--metrics')
        parser.add_argument('--dbg', default=False, action='store_true')

        parser.add_argument('--curio-monitor', default=False, action='store_true',
                            help='Run curio kernel monitor for async debugging.')

        super().__init__(parser, description=desc)

    def postprocess_args(self, args):
        process_graph_args(args)
        self.status = StatusOutput(term=self.term)

    def setup(self, args):
        if args.dbg:
            self.hasher = args.hasher_t(args.ksize)
            self.storage = args.storage.build(*args.storage_args)
            self.streamer_t = args.graph_t
            self.streamer = args.graph_t.build(self.storage, self.hasher)
        else:
            self.streamer_t = libgoetia.StreamHasher[args.hasher_t]
            self.streamer = self.streamer_t.Hasher.build(args.ksize)

        self.ll_processor = self.streamer_t.Processor.build(self.streamer,
                                                            args.interval)

        sample_iter = iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names)
        self.processor = AsyncSequenceProcessor(self.ll_processor, sample_iter)

        self.worker_listener = self.processor.add_listener('worker_q', 'consumer')
        
        if args.metrics:
            self.metrics_stream = AsyncJSONStreamWriter(args.metrics)
            self.worker_listener.on_message(
                Interval,
                StreamHasherRunner.write_metrics_callback,
                self.metrics_stream
            )

        self.worker_listener.on_message(
            Interval,
            lambda msg, status: \
                status.update(msg.t, msg.sequence, msg.seconds_elapsed_interval),
            self.status
        )

        self.worker_listener.on_message(
            SampleStarted,
            lambda msg, status: status.start_sample(msg.sample_name, msg.file_names),
            self.status
        )

        self.worker_listener.on_message(
            SampleFinished,
            lambda msg, status: status.finish_sample(msg.seconds_elapsed_sample),
            self.status
        )

        self.worker_listener.on_message(
            Error,
            lambda msg, status: status.message(
                f'ERROR: {msg.sample_name}->{msg.file_names} at time={msg.t}, sequence={msg.sequence} '\
                f'\n-- BEGIN EXCEPTION --\n{msg.error}\n-- END EXCEPTION --'
            ),
            self.status
        )

    def execute(self, args):
        curio.run(self.processor.start, with_monitor=args.curio_monitor)

    def teardown(self):
        pass

    @staticmethod
    async def write_metrics_callback(msg, stream):
        data = {'t': msg.t,
                'seq_t': msg.sequence,
                'rt_elapsed_interval': msg.seconds_elapsed_interval,
                'rt_elapsed_total': msg.seconds_elapsed_total,
                'sample_name': msg.sample_name}
        await stream.write(data)
