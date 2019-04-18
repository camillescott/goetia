# prometheus.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import socket
from cppyy.gbl import std

def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0


class Instrumentation:

    def __init__(self, port,
                       uri = '/metrics',
                       expose  = True):

        self.Registry = std.make_shared[gbl.prometheus.Registry]()
        ref = std.weak_ptr[gbl.prometheus.Registry](self.Registry.__smartptr__())

        if expose:
            _port = int(port)
            try:
                if is_port_in_use(_port):
                    raise ValueError('Port {0} is in use.'.format(_port))
            except Exception as e:
                raise e

            _address = '127.0.0.1:' + str(port)
            _uri     = uri

            self.Exposer  = std.make_shared[gbl.prometheus.Exposer](_address, _uri);
            self.Exposer.RegisterCollectable(ref)
