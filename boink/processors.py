
import curio
from curio.socket import *

from collections import OrderedDict
import json
import os
import signal
import sys
import threading
from typing import Awaitable, Optional, Callable, Tuple

DEFAULT_SOCKET = '/tmp/boink.sock'

class MessageTypes:

    class Interval:
        @staticmethod
        def generate(name, t, state):
            assert t >= 0
            return 'Interval', {'t': t, 'state': state, 'sample.name': name}
    
    class SampleStarted:
        @staticmethod
        def generate(name):
            return 'SampleStarted', {'sample.name': name}
    
    class SampleFinished:
        @staticmethod
        def generate(name, n):
            return 'SampleFinished', {'sample.name': name, 'N': n}

    class Error:
        @staticmethod
        def generate(name, error):
            return 'Error', {'sample.name': name, 'error': error}

    @staticmethod
    def generate(msg_type, *args):
        type_name, data = msg_type.generate(*args)
        return {'MSG_TYPE': type_name, 'DATA': data}


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
        async for msg in self.q:
            for sub_q in self.subscribers:
                await sub_q.put(msg)
            if msg is None:
                break
        await self.q.task_done()


class AsyncSequenceProcessor:

    def __init__(self, processor,
                       sample_iter,
                       broadcast_socket = DEFAULT_SOCKET):
        """Manages advancing through a concrete FileProcessor
        CRTP subblass asynchronously. The processor pushes Interval
        updates on to the `worker_q`, which are also forwarded
        to an `events_q`. Additional async tasks can subscribe to 
        either queue; the `events_q` is considered the outward-facing
        point.

        `sample_iter` should be conform to that produced by
        `boink.processing.iter_fastx_inputs`.
        
        Args:
            processor (libboink.InserterProcessor<T>): Processor to manage.
            sample_iter (iterator): Iterator over pairs of or single samples.
            broadcast_socket (str, optional): AF_UNIX socket to broadcast
                the events queue on.
        """
        try:
            os.unlink(broadcast_socket)
        except OSError:
            if os.path.exists(broadcast_socket):
                raise
        self.broadcast_socket = broadcast_socket

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

        self.processor = processor
        self.sample_iter = sample_iter
    
    def subscribe(self, channel: str, collection_q: curio.Queue, subscriber_name: str):
        """Subscribe to a queue of the given name to a channel.
        
        Args:
            channel (str): Name of the channel.
            collection_q (curio.Queue): The queue to collect on.
            subscriber_name (str): Name of the subscriber.
        """
        try:
            self.channels[channel].subscribe(collection_q, subscriber_name)
        except KerError:
            print(f'Attmpted to subscribe to invalid channel: "{channel}" does not exist.', file=sys.stderr)
        print(f'{subscriber_name} subscribed to {channel}.', file=sys.stderr)
    
    def unsubscribe(self, channel: str, collection_q: curio.Queue) -> None:
        """Stop receving data from the named channel on the given queue.
        
        Args:
            channel (str): Name of the channel.
            collection_q (curio.Queue): Queue object to remove.
        """
        try:
            self.channels[channel].unsubscribe(collection_q)
        except KerError:
            print(f'Attmpted to unsubscribe from invalid channel: "{channel}" does not exist.')

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
            while True:
                print(f'Starting broadcast server on {self.broadcast_socket}.', file=sys.stderr)
                broadcast_task = await curio.spawn(curio.unix_server,
                                                   self.broadcast_socket,
                                                   self.broadcast_client)
                await restart.get()
                await broadcast_task.cancel()

    @curio.async_thread
    def worker(self) -> None:
        for sample, prefix in self.sample_iter:
            self.worker_q.put(MessageTypes.generate(MessageTypes.SampleStarted, str(sample)))
            try:
                for n_seqs, state in self.processor.chunked_process(*sample):
                    self.worker_q.put(MessageTypes.generate(MessageTypes.Interval, str(sample), n_seqs, state.get()))
            except Exception as e:
                self.worker_q.put(MessageTypes.generate(MessageTypes.Error, str(sample), str(e)))
                raise
            self.worker_q.put(MessageTypes.generate(MessageTypes.SampleFinished, str(sample), n_seqs))

    async def status(self):
        """Writes everything from the events queue to `sys.stderr`.
        """
        try:
            msg_q = curio.Queue()
            self.subscribe('events_q', msg_q, 'status')

            while True:
                msg = await msg_q.get()

                if msg is None:
                    await msg_q.task_done()
                    break

                print(json.dumps(msg), file=sys.stderr)
                await msg_q.task_done()

        except curio.CancelledError:
            raise
        finally:
            self.unsubscribe('events_q', msg_q)

    async def run(self, extra_tasks = None) -> None:
        async with curio.TaskGroup() as g:

            # Each channel has its own dispatch task
            # to send data to its subscribers
            for channel_name, channel in self.channels.items():
                await g.spawn(channel.dispatch)

            #await g.spawn(self.broadcaster)
            status = await g.spawn(self.status)

            # Spawn extra tasks to run
            if extra_tasks is not None:
                for task in extra_tasks:
                    await g.spawn(task)

            # and now we spawn the worker to iterate through
            # the processor and wait for it to finish
            await self.worker()
            await g.cancel_remaining()