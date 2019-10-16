#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : pythonize_signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.10.2019

from numpy import zeros
import numpy

import json
import typing

from boink.pythonizors import utils

def pythonize_boink_signatures(klass, name):
    if name == 'SourmashSignature':

        def to_sourmash(self):
            try:
                from sourmash import MinHash
            except ImportError:
                print('Must install sourmash to conver to sourmash.MinHash',
                      file=sys.stderr)
                return None

            sig = MinHash(self.num, self.ksize)
            sig.add_many(self.mins)

            return sig

        klass.Signature.to_sourmash = to_sourmash

    is_inst, template =  utils.is_template_inst(name, 'UnikmerSignature')
    if is_inst:
        def signature(self) -> numpy.ndarray:
            """signature

            Returns:
                numpy.ndarray: Numpy array with the signature vector.
            """
            sig = zeros(len(self))
            raw = self.get_signature()
            for i, count in enumerate(raw):
                sig[i] = count
            return sig

        def __len__(self) -> int:
            return self.get_size()

        def to_dict(self, name: str) -> dict:
            """to_dict

            Args:
                name (str): Convert the signature metadata and signature to a dictionary.

            Returns:
                dict: Signature metadata.
            """
            data = {'signature': self.signature.tolist(),
                    'W'        : self.K,
                    'K'        : self.bucket_K,
                    'size'     : len(self),
                    'name'     : name}
            return data

        def save(self, stream: typing.TextIO, name: str) -> None:
            """Save the signature to disk.

            Args:
                stream (file): File handle to save to.
                name (str): Name of the signature.
            """
            data = [self.to_dict(name)]

            json.dump(data, stream)

        klass.Signature.signature = property(signature)
        klass.Signature.__len__   = __len__
        klass.Signature.to_dict   = to_dict
        klass.Signature.save      = save
