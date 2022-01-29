#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : processors.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 05.11.2021

from collections import OrderedDict, defaultdict
from enum import Enum, unique as unique_enum
import functools
import inspect
import json
import signal
import sys
from typing import Any

import curio

from goetia import libgoetia
from goetia.messages import (
    AllMessages,
    EndStream,
    Error,
    Interval,
    SampleFinished,
    SampleStarted,
)
from goetia.utils import is_iterable


DEFAULT_SOCKET = '/tmp/goetia.sock'
DEFAULT_INTERVAL = libgoetia.IntervalCounter.DEFAULT_INTERVAL


class QueueManager:

    def __init__(self, q: curio.UniversalQueue, name: str):
        assert type(q) is curio.UniversalQueue
        self.q = q
        self.name = name
        self.subscribers = set()
        self.subscriber_names = {}
    
    def subscribe(self, q: curio.UniversalQueue, name: str):
        if q not in self.subscribers:
            self.subscribers.add(q)
            self.subscriber_names[q] = name

    def unsubscribe(self, q: curio.UniversalQueue) -> None:
        try:
            self.subscribers.remove(q)
            del self.subscriber_names[q]
        except:
            pass

    async def kill(self) -> None:
        await self.q.put(None)
    
    async def dispatch(self) -> None:
        while True:
            msg = await self.q.get()
            for sub_q in self.subscribers:
                await sub_q.put(msg)
            await self.q.task_done()

            if msg is None:
                break


class MessageHandler:
    
    def __init__(self, name, subscription):
        self.subscription = subscription
        self.name = name
        self.handlers = defaultdict(list)
    
    async def task(self, error_handler = None):
        try:
            msg_q = curio.Queue()
            self.subscription.subscribe(msg_q, self.name)
            while True:
                msg = await msg_q.get()
                if msg is None:
                    await msg_q.task_done()
                    break
                
                for callback, args, kwargs in self.handlers[type(msg)]:
                    if inspect.iscoroutinefunction(callback):
                        await callback(msg, *args, **kwargs)
                    else:
                        callback(msg, *args, **kwargs)

                for callback, args, kwargs in self.handlers[AllMessages]:
                    if inspect.iscoroutinefunction(callback):
                        await callback(msg, *args, **kwargs)
                    else:
                        callback(msg, *args, **kwargs)

                await msg_q.task_done()
        except curio.CancelledError:
            raise
        except Exception as e:
            if error_handler is not None:
                error_handler(e)
            else:
                raise
        else:
            self.subscription.unsubscribe(msg_q)
    
    def on_message(self, msg_class, callback, *args, **kwargs):
        assert type(msg_class) is type
        self.handlers[msg_class].append((callback, args, kwargs))


@unique_enum
class RunState(Enum):
    READY = 0
    RUNNING = 1
    SIGINT = 2
    STOP_SATURATED = 3
    STOP_ERROR = 4
    STOP = 5


class AsyncSequenceProcessor:

    def __init__(self, processor,
                       sample_iter,
                       echo = None,
                       broadcast_socket = None):
        """Manages advancing through a concrete FileProcessor
        subblass asynchronously. The processor pushes Interval
        updates on to the `worker_q`, which are also forwarded
        to an `events_q`. Additional async tasks can subscribe to 
        either queue; the `events_q` is considered the outward-facing
        point.

        `sample_iter` should be conform to that produced by
        `goetia.processing.iter_fastx_inputs`.
        
        Args:
            processor (libgoetia.InserterProcessor<T>): Processor to manage.
            sample_iter (iterator): Iterator over pairs of or single samples.
            echo (bool): Whether to echo `events_q` to the terminal.
            broadcast_socket (str, optional): AF_UNIX socket to broadcast
                the events queue on.
        """
 
        self.worker_q = curio.UniversalQueue()
        self.worker_subs = QueueManager(self.worker_q, 'worker_q')

        self.events_q = curio.UniversalQueue()
        self.events_subs = QueueManager(self.events_q, 'events_q')

        self.channels = OrderedDict()
        self.channels[self.worker_subs.name] = self.worker_subs
        self.channels[self.events_subs.name] = self.events_subs

        # We want everything from the worker q to also end
        # up on the events q
        self.subscribe('worker_q', self.events_q, 'events_q')

        self.listener_tasks = []

        self.processor = processor
        self.sample_iter = sample_iter

        self.run_echo = echo is not None
        self.echo_file = '/dev/stderr' if echo is True else echo

        self.state = RunState.READY
        self.processed = set()

        #super().__init__(broadcast_socket)
    
    def get_channel(self, channel: str) -> QueueManager:
        """Query for the given channel name.
        
        Args:
            channel (str): The channel name.
        
        Returns:
            QueueManager: Manager for the channel.
        """        
        try:
            return self.channels[channel]
        except KeyError:
            print(f'Requested invalid channel: "{channel}" does not exist.', file=sys.stderr)
            raise
    
    def subscribe(self, channel_name: str, 
                        collection_q: curio.UniversalQueue,
                        subscriber_name: str) -> None:
        """Subscribe a queue of the given name to a channel.
        
        Args:
            channel_name (str): Name of the channel.
            collection_q (curio.Queue): The queue to collect on.
            subscriber_name (str): Name of the subscriber.
        """
        self.get_channel(channel_name).subscribe(collection_q, subscriber_name)
    
    def unsubscribe(self, channel_name: str,
                          collection_q: curio.UniversalQueue) -> None:
        """Stop receving data from the named channel on the given queue.
        
        Args:
            channel_name (str): Name of the channel.
            collection_q (curio.Queue): Queue object to remove.
        """
        self.get_channel(channel_name).unsubscribe(collection_q)
        
    def add_listener(self, channel_name: str,
                           subscriber_name: str) -> MessageHandler:
        channel = self.get_channel(channel_name)
        listener = MessageHandler(subscriber_name, channel)
        self.listener_tasks.append(listener.task)
        return listener

    def worker(self) -> None:
        time, n_seqs = 0, 0
        for sample, name in self.sample_iter:
            self.worker_q.put(SampleStarted(sample_name=name,  # type: ignore
                                            file_names=sample,
                                            t=time,
                                            sequence=n_seqs))
            try:
                for n_seqs, time, n_skipped in self.processor.chunked_process(*sample):
                    
                    if self.state is RunState.STOP_SATURATED:
                        # Saturation is tripped externally: just return immediately.
                        return

                    if self.state is RunState.SIGINT:
                        # If we're interrupted, inform our listeners that something went wrong.
                        self.worker_q.put(Error(t=time,  # type: ignore
                                                sequence=n_seqs,
                                                sample_name=name,
                                                file_names=sample,
                                                error='Process terminated (SIGINT).'))
                        return

                    self.worker_q.put(Interval(t=time,  # type: ignore
                                               sequence=n_seqs,
                                               sample_name=name, 
                                               file_names=sample))

                self.processed.add(tuple(sample))
                self.worker_q.put(SampleFinished(t=time,  # type: ignore
                                                 sequence=n_seqs,
                                                 sample_name=name,
                                                 file_names=sample))
            except Exception as e:
                self.worker_q.put(Error(t=time,  # type: ignore
                                        sequence=n_seqs,
                                        sample_name=name,
                                        file_names=sample,
                                        error=str(e.__traceback__))) 
                return
            finally:
                self.worker_q.put(EndStream(t=time,  # type: ignore
                                            sequence=n_seqs))
    
    #def on_error(self, exception):
    #    self.worker_q.put(Error(t=self.processor.time_elapsed(),
    #                            sequence=n_seqs,
    #                            sample_name=name,
    #                            error=f'At sequence {self.processor.n_sequences()}: {str(e)}',
    #                            file_names=sample))
    #    self.state = RunState.STOP_ERROR

    async def start(self, extra_tasks = None) -> None:
        try:
            async with curio.TaskGroup() as g:

                # each channel has its own dispatch task
                # to send data to its subscribers
                for channel_name, channel in self.channels.items():
                    await g.spawn(channel.dispatch)

                # start up AF_UNIX broadcaster if desired
                # if self.broadcast_socket is not None:
                #    await g.spawn(self.broadcaster)

                if self.run_echo:
                    listener = self.add_listener('events_q', 'echo')
                    async def echo(msg):
                        mode = 'w' if self.echo_file in ['/dev/stdout', '/dev/stderr'] else 'a'
                        async with curio.aopen(self.echo_file, mode) as fp:
                            await fp.write(f'{msg.to_yaml()}\n')

                    listener.on_message(AllMessages,
                                        echo)
                
                # spawn tasks from listener callbacks
                for task in self.listener_tasks:
                    await g.spawn(task)

                # spawn extra tasks to run
                if extra_tasks is not None:
                    for task in extra_tasks:
                        await g.spawn(task)

                # and now we spawn the worker to iterate through
                # the processor and wait for it to finish
                self.state = RunState.RUNNING
                signal.signal(signal.SIGINT, lambda signo, frame: self.interrupt())

                # give just a bit of time for the listeners to all spin up
                await curio.sleep(0.05)
                # then spawn the worker
                w = await g.spawn_thread(self.worker)
                await w.join()
                await curio.sleep(0.05)

                await self.worker_subs.kill()
        except Exception as e:
            print(e, file=sys.stderr)

    def stop(self) -> None:
        self.state = RunState.STOP

    def interrupt(self) -> None:
        self.state = RunState.SIGINT

    def saturate(self) -> None:
        self.state = RunState.STOP_SATURATED


def every_n_intervals(func, n=1):
    poller = libgoetia.metrics.IntervalCounter(n)
    @functools.wraps(func)
    async def wrapped(msg, *args, **kwargs):
        assert isinstance(msg, Interval)
        if poller.poll():
            await func(msg, *args, **kwargs)
    return wrapped


class AsyncJSONStreamWriter:

    def __init__(self, filename: str):
        '''Asynchronously write JSON data to a file.

        Writes a stream of JSON objects to a file. The top-level
        element is always a list; list items can be any valid JSON
        type.

        Args:
            filename: Path of the target file to write to.
        '''
        self.filename = filename
        self.n_writes = 0
        
        with open(self.filename, 'w') as fp:
            fp.write('[')

    def __del__(self):
        with open(self.filename, 'a') as fp:
            fp.write(']')

    async def write(self, data: Any, expand: bool = True):
        '''Write the given data is a JSON element to the stream.
        Strings will be written assuming they are already valid JSON;
        this could result in malformed JSON, so care must be taken.
        Other data types are passed to json.dumps for serialization.

        Args:
            data: Data to coerce to JSON.
            expand: If True, iterables will be expanded into the stream
                    rather than appended as a single item.
        '''

        buf = ''

        async with curio.aopen(self.filename, 'a') as fp:
            if self.n_writes != 0:
                await fp.write(',\n')
            if isinstance(data, str):
                # assume already valid JSON object
                buf = data
            elif expand and is_iterable(data) and not isinstance(data, dict):
                # extend the top level list rather than
                # adding the iterable as an item
                buf = ','.join((json.dumps(item) for item in data))
            else:
                buf = json.dumps(data)
            await fp.write(buf)

        self.n_writes += 1

