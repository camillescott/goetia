
import curio
from curio.socket import *

from collections import OrderedDict
import json
import os
import signal
import sys
import threading
from typing import Awaitable, Optional, Callable, Tuple

from boink.messages import *

DEFAULT_SOCKET = '/tmp/boink.sock'


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


class MessageHandler:
    
    def __init__(self, name, subscription):
        self.subscription = subscription
        self.name = name
        self.handlers = {}
    
    async def task(self):
        try:
            msg_q = curio.Queue()
            self.subscription.subscribe(msg_q, self.name)
            while True:
                msg = await msg_q.get()
                if msg is None:
                    await msg_q.task_done()
                    break
                try:
                    callback, args = self.handlers[type(msg)]
                except KeyError:
                    pass
                else:
                    callback(msg, *args)

                if AllMessages in self.handlers:
                    callback, args = self.handlers[AllMessages]
                    callback(msg, *args)

                await msg_q.task_done()
        except curio.CancelledError:
            raise
        finally:
            self.subscription.unsubscribe(msg_q)
    
    def on_message(self, msg_class, callback, *args):
        assert type(msg_class) is type
        self.handlers[msg_class] = (callback, args)


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


class AsyncSequenceProcessor(UnixBroadcasterMixin):

    def __init__(self, processor,
                       sample_iter,
                       echo = True,
                       broadcast_socket = None):
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
        self.run_echo = echo

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
        """Subscribe to a queue of the given name to a channel.
        
        Args:
            channel_name (str): Name of the channel.
            collection_q (curio.Queue): The queue to collect on.
            subscriber_name (str): Name of the subscriber.
        """
        self.get_channel(channel_name).subscribe(collection_q, subscriber_name)
        print(f'{subscriber_name} subscribed to {channel_name}.', file=sys.stderr)
    
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

    @curio.async_thread
    def worker(self) -> None:
        for sample, prefix in self.sample_iter:
            self.worker_q.put(SampleStarted(sample_name=str(sample)))
            try:
                for n_seqs, state in self.processor.chunked_process(*sample):
                    if not state.end:
                        self.worker_q.put(Interval(t=n_seqs,
                                                   sample_name=str(sample), 
                                                   state=state.get()))
                self.worker_q.put(SampleFinished(t=n_seqs,
                                                 sample_name=str(sample)))
            except Exception as e:
                self.worker_q.put(Error(t=n_seqs,
                                        sample_name=str(sample),
                                        error=str(e)))
                raise

    async def run(self, extra_tasks = None) -> None:
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
                    listener.on_message(AllMessages,
                                        lambda msg: print(msg.to_json(), file=sys.stderr))
                
                # spawn tasks from listener callbacks
                for task in self.listener_tasks:
                    await g.spawn(task)

                # spawn extra tasks to run
                if extra_tasks is not None:
                    for task in extra_tasks:
                        await g.spawn(task)

                # and now we spawn the worker to iterate through
                # the processor and wait for it to finish
                await self.worker()
                await self.worker_subs.kill()
                #await g.cancel_remaining()
        except curio.TaskGroupError as e:
            print(f'TaskGroup failed: {e.errors}')
            for task in e:
                print(f'Task {task} failed due to: {task.exception}')