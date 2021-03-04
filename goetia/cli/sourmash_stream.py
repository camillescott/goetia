#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : sourmash_stream.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 04.06.2020

from sourmash import SourmashSignature, save_signatures

from goetia.sketches import SourmashSketch
from goetia.cli.cli import format_filenames
from goetia.cli.signature_runner import SignatureRunner


class SourmashRunner(SignatureRunner):

    def __init__(self, parser):

        parser.add_argument('-K', default=31, type=int)
        parser.add_argument('-N', default=2500, type=int)
        parser.add_argument('--scaled', default=0, type=int)

        super().__init__(parser)

    @staticmethod
    def _make_signature(args):
        if args.scaled:
            return SourmashSketch.Signature.build(0, args.K, False, False, False, 42, args.scaled)
        else:
            return SourmashSketch.Signature.build(args.N, args.K, False, False, False, 42, 0)

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
    def _save_signatures(sigs, args):
        with open(args.save_sig, 'w') as fp:
            save_signatures(sigs, fp=fp)
