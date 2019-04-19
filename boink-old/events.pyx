# events.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref

cdef class EventListener:

    def __cinit__(self, *args, **kwargs):
        pass

    @staticmethod
    cdef EventListener _wrap(shared_ptr[_EventListener] listener):
        cdef EventListener el = EventListener()
        el._listener = listener
        return el

    def stop(self):
        deref(self._listener).exit_thread()

    def wait_on_processing(self, int min_events=0):
        deref(self._listener).wait_on_processing(min_events)


cdef class EventNotifier:

    def __cinit__(self, *args, **kwargs):
        pass

    @staticmethod
    cdef EventNotifier _wrap(shared_ptr[_EventNotifier] notifier):
        cdef EventNotifier en = EventNotifier()
        en._notifier = notifier 
        return en

    def register_listener(self, EventListener listener):
        deref(self._notifier).register_listener(listener._listener.get())

    def stop_listeners(self):
        deref(self._notifier).stop_listeners()

    def clear_listeners(self):
        deref(self._notifier).clear_listeners()
