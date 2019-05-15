# prometheus.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import socket
from cppyy import gbl
from cppyy.gbl import std


def get_prometheus_args(parser):
    group = parser.add_argument_group('prometheus')
    group.add_argument('--port', default=None,
                        help='Port to expose prometheus metrics.')
    return group


def print_prometheus_args(args):
    print('* Exposing prometheus metrics on port', args.port, file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)


def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0


class Instrumentation:

    def __init__(self, port = None,
                       uri = '/metrics',
                       expose  = False):

        self.Registry = std.make_shared[gbl.prometheus.Registry]()
        ref = std.weak_ptr[gbl.prometheus.Registry](self.Registry.__smartptr__())

        if expose:
            if port is None:
                raise ValueError('Must specificy prometheus port to expose.')
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
