#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : pythonize_signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.10.2019

import sys

import cppyy.ll
import numpy as np

from goetia.pythonizors import utils


def pythonize_goetia_signatures(klass, name):
    if name == 'SourmashSignature':

        def to_sourmash(self):
            try:
                from sourmash import MinHash
            except ImportError:
                print('Must install sourmash to convert to sourmash.MinHash',
                      file=sys.stderr)
                return None

            sig = MinHash(self.num(),
                          self.ksize(),
                          is_protein=self.is_protein(),
                          dayhoff=self.dayhoff(),
                          hp=self.hp(),
                          track_abundance=self.track_abundance(),
                          seed=self.seed(),
                          mins=self.mins(),
                          max_hash=self.max_hash())

            return sig

        klass.Signature.to_sourmash = to_sourmash

    is_inst, template =  utils.is_template_inst(name, 'UnikmerSignature')
    if is_inst:
        def to_numpy(self) -> np.ndarray:
            """signature

            Returns:
                numpy.ndarray: Numpy array with the signature vector.
            """
            buffer = self.get_sketch_as_buffer()
            buffer.reshape((len(self),))
            return np.frombuffer(buffer, dtype=np.uint64, count=len(self))
        
        def __len__(self) -> int:
            return self.get_size()

        klass.Signature.to_numpy = to_numpy
        klass.Signature.__len__   = __len__

        def wrap_build(build_func):
            def wrapped(W, K, ukhs=None, storage_args=None):
                if ukhs is None:
                    ukhs = klass.ukhs_type.load(W, K)
                if storage_args is None:
                    sig = build_func(W, K, ukhs.__smartptr__())
                else:
                    sig = build_func(W, K, ukhs.__smartptr__(), klass.storage_type.make_params(*storage_args))
                return sig
            return wrapped
        
        klass.Signature.build = wrap_build(klass.Signature.build)
