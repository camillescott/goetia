#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : sourmash_stream.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 04.06.2020

from sourmash import SourmashSignature, save_signatures
from sourmash.utils import rustcall, decode_str
from sourmash._lowlevel import ffi, lib

from goetia.sketches import SourmashSketch
from goetia.cli.cli import format_filenames
from goetia.cli.signature_runner import SignatureRunner


desc = '''
{term.italic}sketch module: {term.normal}sourmash stream

    Build a sourmash minhash sketch in a streaming manner.
    Between-signature distances normalized by number of 
    observed {term.italic}k{term.normal}-mers are reported during construction,
    and optionally, a saturation metric can be supplied to cut the
    construction off early.
'''


class SourmashRunner(SignatureRunner):

    def __init__(self, parser):

        parser.add_argument('-K', default=31, type=int)
        parser.add_argument('-N', default=2500, type=int)
        parser.add_argument('--scaled', default=0, type=int)

        super().__init__(parser, description=desc)

    @staticmethod
    def _make_signature(args):
        if args.scaled:
            return SourmashSketch.Sketch.build(0, args.K, False, False, False, 42, args.scaled)
        else:
            return SourmashSketch.Sketch.build(args.N, args.K, False, False, False, 42, 0)

    @staticmethod
    def _make_processor(signature, args):
        return SourmashSketch.Processor.build(signature,
                                              args.interval)

    def _distance_func(self, sigs):
        sig_a, sig_b = sigs
        sim = sig_a.similarity(sig_b)
        return sim
    
    @staticmethod
    def _convert_signature(sig, msg):
        return SourmashSignature(sig.to_sourmash(),
                                 name=f'{msg.sample_name}:{msg.t}',
                                 filename=format_filenames(msg.file_names))        
    
    @staticmethod
    def _save_signature(sig, args):
        with open(args.save_sig, 'w') as fp:
            save_signatures([sig], fp=fp)

    @staticmethod
    def _serialize_signature(sig, args):
        return decode_str(rustcall(lib.signature_save_json, sig._objptr))
