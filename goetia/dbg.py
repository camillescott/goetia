#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : dbg.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

from goetia import libgoetia
from goetia.hashing import typenames as hasher_types
from goetia.storage import get_storage_args, process_storage_args


dBG = libgoetia.dBG
STATES = libgoetia.TraversalState

def Graph(storage_type, shifter_type):

    def build(K, storage_args, shifter_args):
        storage = storage_type.build(*storage_args)
        shifter = shifter_type.build(K, *shifter_args)
        graph_type = dBG[storage_type, shifter_type]

        return graph_type.build(storage, hasher)

    return build


def get_graph_args(parser):
    group = parser.add_argument_group('dBG')
    group.add_argument('-K', '--ksize', type=int, default=31)
    group.add_argument('-H', '--hasher', choices=[name for t, name in hasher_types],
                       default='FwdLemireShifter')

    get_storage_args(parser)

    return group


def process_graph_args(args):
    process_storage_args(args)

    args.hasher_t = getattr(libgoetia.hashing, args.hasher)
    args.graph_t = libgoetia.dBG[args.storage, args.hasher_t]


def print_dBG_args(args):
    print('* dBG will be order', args.ksize, file=sys.stderr)
    print('* dBG will have underlying storage of', args.storage, file=sys.stderr)
    print('*','*' * 10, '*', sep='\n', file=sys.stderr)
