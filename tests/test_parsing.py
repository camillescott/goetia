# goetia/tests/test_parsing.py
# Copyright (C) 2020 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from .utils import *

from goetia.parsing import FastxParser, SplitPairedReader
from goetia.alphabets import DNA_SIMPLE, DNAN_SIMPLE, IUPAC_NUCL

alphabets = [DNA_SIMPLE, DNAN_SIMPLE, IUPAC_NUCL]

def test_parser_sequence(random_fasta):
    sequences, path = random_fasta(10)
    print(path)
    parser = FastxParser[DNA_SIMPLE].build(path)

    parsed = [record.sequence for record in parser]
    assert parsed == sequences


def test_parser_name(random_fasta):
    sequences, path = random_fasta(10)
    parser = FastxParser[DNA_SIMPLE].build(path)

    parsed = [record.name for record in parser]
    assert parsed == [str(i) for i in range(10)]


@pytest.mark.parametrize('alphabet', alphabets)
def test_validation_uppercase(fastx_writer, alphabet):
    sequences = ['aaacctggtt', 'cctggagggg']
    path = fastx_writer(sequences)
    parser = FastxParser[alphabet].build(str(path))

    parsed = [record.sequence for record in parser]
    assert parsed == [sequence.upper() for sequence in sequences]


def test_validation_bad_char_no_strict(fastx_writer):
    sequences = ['AAAAAANA']
    path = fastx_writer(sequences)
    parser = FastxParser[DNA_SIMPLE].build(str(path))
    parsed = list(parser)

    assert len(parsed) == 0
    assert parser.n_skipped() == 1


@pytest.mark.parametrize('alphabet', alphabets)
def test_empty_file(alphabet, fastx_writer):
    path = fastx_writer([])
    parser = FastxParser[alphabet].build(str(path))

    assert [record.sequence for record in parser] == []


@pytest.mark.parametrize('alphabet', alphabets)
def test_empty_file(alphabet):
    with pytest.raises(Exception):
        parser = FastxParser[alphabet].build('bad_path.fa')
        list(parser)


def test_split_paired_read(datadir):
    left = datadir('left.fq')
    right = datadir('right.fq')

    seqs = []
    for left, right in SplitPairedReader[FastxParser[DNA_SIMPLE]].build(left, right):
        seqs.append((left, right))

    print(seqs)
    assert len(seqs) == 25
