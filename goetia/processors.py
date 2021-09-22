
import curio
from curio.socket import *

from collections import defaultdict, OrderedDict
from enum import Enum, unique as unique_enum
import functools
import inspect
import json
import os
import signal
import sys
import threading
from typing import Awaitable, Optional, Callable, Tuple

from goetia.messages import (Interval, DistanceCalc, SampleStarted, SampleFinished,
                             SampleSaturated, Error, AllMessages)
from goetia import libgoetia


DEFAULT_SOCKET = '/tmp/goetia.sock'
DEFAULT_INTERVAL = libgoetia.metrics.IntervalCounter.DEFAULT_INTERVAL


class QueueManager:

    def __init__(self, q: curio.UniversalQueue, name: str):
        assert type(q) is curio.UniversalQueue
        self.q = q
        self.name = name
        self.subscribers = set()
        self.subscriber_names = {}
    
    def subscribe(self, q: curio.Queue, name: str):
        if q not in self.subscribers:
            self.subscribers.add(q)
            self.subscriber_names[q] = name

    def unsubscribe(self, q: curio.Queue) -> None:
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
        finally:
            self.subscription.unsubscribe(msg_q)
    
    def on_message(self, msg_class, callback, *args, **kwargs):
        assert type(msg_class) is type
        self.handlers[msg_class].append((callback, args, kwargs))


class UnixBroadcasterMixin:

    def __init__(self, broadcast_socket = DEFAULT_SOCKET):
        if broadcast_socket is not None:
            try:
                os.unlink(broadcast_socket)
            except OSError:
                if os.path.exists(broadcast_socket):
                    raise
        self.broadcast_socket = broadcast_socket

    async def broadcast_client(self, client: curio.io.Socket, addr: Tuple[str, int]) -> None:
        client_name = hash(client) # i guess getpeername() doesn't work with AF_UNIX
        print(f'Unix socket connection: {client_name}', file=sys.stderr)

        stream = client.as_stream()
        bcast_q = curio.Queue()
        self.subscribe('events_q', bcast_q, f'broadcast_client:{client_name}')
        
        try:
            while True:
                block = await bcast_q.get()

                string = json.dumps(block) + '\n'
                await curio.timeout_after(60, stream.write, string.encode('ascii'))
        except curio.CancelledError:
            await stream.write(json.dumps([(0, 'END_STREAM', -1)]).encode('ascii'))
            raise
        except (BrokenPipeError, curio.TaskTimeout):
            print(f'Unix socket closed: {client_name}', file=sys.stderr)
        finally:
            self.unsubscribe('events_q', bcast_q)

    async def broadcaster(self) -> None:
        async with curio.SignalQueue(signal.SIGHUP) as restart:
            if self.broadcast_socket is not None:
                while True:
                    print(f'Starting broadcast server on {self.broadcast_socket}.', file=sys.stderr)
                    broadcast_task = await curio.spawn(curio.unix_server,
                                                    self.broadcast_socket,
                                                    self.broadcast_client)
                    await restart.get()
                    await broadcast_task.cancel()


@unique_enum
class RunState(Enum):
    READY = 0
    RUNNING = 1
    SIGINT = 2
    STOP_SATURATED = 3
    STOP_ERROR = 4
    STOP = 5


class AsyncSequenceProcessor(UnixBroadcasterMixin):

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

        super().__init__(broadcast_socket)
    
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
                        collection_q: curio.Queue,
                        subscriber_name: str) -> None:
        """Subscribe a queue of the given name to a channel.
        
        Args:
            channel_name (str): Name of the channel.
            collection_q (curio.Queue): The queue to collect on.
            subscriber_name (str): Name of the subscriber.
        """
        self.get_channel(channel_name).subscribe(collection_q, subscriber_name)
    
    def unsubscribe(self, channel_name: str,
                          collection_q: curio.Queue) -> None:
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
        for sample, name in self.sample_iter:
            self.worker_q.put(SampleStarted(sample_name=name, file_names=sample))
            try:
                for n_seqs, time, n_skipped in self.processor.chunked_process(*sample):
                    
                    if self.state is RunState.STOP_SATURATED:
                        # Saturation is tripped externally: just return immediately.
                        return

                    if self.state is RunState.SIGINT:
                        # If we're interrupted, inform our listeners that something went wrong.
                        self.worker_q.put(Error(t=time,
                                                sequence=n_seqs,
                                                sample_name=name,
                                                file_names=sample,
                                                error='Process terminated (SIGINT).'))
                        return

                    self.worker_q.put(Interval(t=time,
                                               sequence=n_seqs,
                                               sample_name=name, 
                                               file_names=sample))

                self.processed.add(tuple(sample))
                self.worker_q.put(SampleFinished(t=time,
                                                 sequence=n_seqs,
                                                 sample_name=name,
                                                 file_names=sample))
            except Exception as e:
                self.worker_q.put(Error(t=time,
                                        sequence=n_seqs,
                                        sample_name=name,
                                        file_names=sample,
                                        error=e.__traceback__))
                return
    
    def on_error(self, exception):
        self.worker_q.put(Error(t=self.processor.time_elapsed(),
                                sequence=n_seqs,
                                sample_name=name,
                                error=f'At sequence {self.processor.n_sequences()}: {str(e)}',
                                file_names=sample))
        self.state = RunState.STOP_ERROR

    async def start(self, extra_tasks = None) -> None:
        try:
            async with curio.TaskGroup() as g:

                # each channel has its own dispatch task
                # to send data to its subscribers
                for channel_name, channel in self.channels.items():
                    await g.spawn(channel.dispatch)

                # start up AF_UNIX broadcaster if desired
                if self.broadcast_socket is not None:
                    await g.spawn(self.broadcaster)

                if self.run_echo:
                    listener = self.add_listener('events_q', 'echo')
                    async def echo(msg):
                        mode = 'w' if self.echo_file in ['/dev/stdout', '/dev/stderr'] else 'a'
                        async with curio.aopen(self.echo_file, mode) as fp:
                            await fp.write(f'{msg.to_json()}\n')

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


def at_modulo_interval(func, modulus=1):
    @functools.wraps(func)
    async def wrapped(msg, *args, **kwargs):
        assert isinstance(msg, Interval)
        if msg.t % modulus == 0:
            await func(msg, *args, **kwargs)
    return wrapped
