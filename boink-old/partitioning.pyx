# boink/partitioning.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from __future__ import print_function
import argparse
import itertools
import json
import os
import sys

import cython
from cython.operator cimport dereference as deref, preincrement as inc

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr, make_shared

from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.graphs cimport Countgraph, Nodegraph, Hashgraph
from khmer.khmer_args import build_counting_args, create_countgraph
from khmer.khmer_logger import (configure_logging, log_info, log_error,
                                log_warn)

from khmer._oxli.partitioning cimport (CpStreamingPartitioner,
                                       StreamingPartitioner, Component,
                                       ComponentPtr)
from khmer._oxli.parsing cimport BrokenPairedReader, SplitPairedReader, FastxParser
from khmer._oxli.parsing import BrokenPairedReader, SplitPairedReader, FastxParser
from khmer._oxli.sequence cimport Sequence
from khmer._oxli.sequence import Sequence
from khmer._oxli.oxli_types cimport *

from boink.stats cimport PartitionFunction, PartitionCoverage, PartitionCoverageSlice
from boink.utils cimport _bstring

cdef class ConditionalPartitioner(StreamingPartitioner):

    def __cinit__(self, Hashgraph graph, tag_density=None, info_filename=None):
        self.min_fitness = 0.95
        self.info_filename = info_filename
        if self.info_filename is not None:
            self.info_fp = open(self.info_filename, 'w')
            self.info_fp.write('read_n,fitness,n_tags,n_merged,comp_size,comp_cov,root_id\n')
        else:
            self.info_fp = None

    def __dealloc__(self):
        try:
            self.info_fp.close()
        except Exception:
            pass

    def consume(self, Sequence sequence, PartitionFunction func=None):
        if func is None:
            return deref(self._this).consume(sequence._obj.sequence)
        cdef float fitness
        self._consume_conditional(sequence._obj.sequence, func, fitness)
        return fitness

    def consume_pair(self, Sequence first, Sequence second,
                       PartitionFunction func=None):

        if func is None:
            return deref(self._this).consume_pair(first._obj.sequence,
                                                    second._obj.sequence)
        cdef float fitness
        self._consume_conditional_pair(first._obj.sequence, second._obj.sequence,
                                         func, fitness)
        return fitness
        
    cdef int _consume_conditional(self, string sequence, PartitionFunction func, 
                                  float& fitness) except -1:

        cdef vector[HashIntoType] tags
        cdef KmerQueue seeds
        cdef set[HashIntoType] seen

        cdef uint64_t n_new = deref(self._this).seed_sequence(sequence, tags, seeds, seen)
        self.n_consumed += 1
        (&fitness)[0] = func._evaluate_tags(sequence, tags)
        if fitness < self.min_fitness: # blerp
            if self.info_fp is not None:
                self.info_fp.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(self.n_consumed,
                                                                          fitness,
                                                                          tags.size(),
                                                                          None,
                                                                          None,
                                                                          None,
                                                                          None))
            return n_new

        deref(self._this).find_connected_tags(seeds, tags, seen, False)
        cdef uint32_t n_merged = deref(self._this).create_and_merge_components(tags)
        cdef ComponentPtr compptr
        cdef float comp_cov = -1
        compptr = deref(self._this).get(deref(tags.begin()))
        if n_merged > 1:
            comp_cov = Component._mean_tag_count(compptr,
                                                 self.graph._hg_this.get())
        if self.info_fp is not None:
            self.info_fp.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(self.n_consumed,
                                                                      fitness,
                                                                      tags.size(),
                                                                      n_merged,
                                                                      deref(compptr).get_n_tags(),
                                                                      comp_cov,
                                                                      deref(compptr).component_id))
        return n_new

    cdef int _consume_conditional_pair(self, string first, string second,
                                         PartitionFunction func, float& fitness) except -1:

        cdef vector[HashIntoType] tags
        cdef KmerQueue seeds
        cdef set[HashIntoType] seen
        cdef uint64_t n_new = deref(self._this).seed_sequence(first, tags, seeds, seen)
        n_new +=  deref(self._this).seed_sequence(second, tags, seeds, seen)
        self.n_consumed += 2

        (&fitness)[0] = func._evaluate_tags(first, tags)
        if fitness < self.min_fitness: # blerp
            if self.info_fp is not None:
                self.info_fp.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(self.n_consumed,
                                                                          fitness,
                                                                          tags.size(),
                                                                          None,
                                                                          None,
                                                                          None,
                                                                          None))
            return n_new
        
        deref(self._this).find_connected_tags(seeds, tags, seen, False)
        cdef uint32_t n_merged = deref(self._this).create_and_merge_components(tags)
        cdef ComponentPtr compptr
        cdef float comp_cov = -1
        compptr = deref(self._this).get(deref(tags.begin()))
        if n_merged > 1:
            comp_cov = Component._mean_tag_count(compptr,
                    self.graph._hg_this.get())
        if self.info_fp is not None:
            self.info_fp.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(self.n_consumed,
                                                                      fitness,
                                                                      tags.size(),
                                                                      n_merged,
                                                                      deref(compptr).get_n_tags(),
                                                                      comp_cov,
                                                                      deref(compptr).component_id))
        return n_new


def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


cdef class DoConditionalPartitioning:

    def __init__(self, args=sys.argv[1:]):
        self.args = self.parse_args(args)
        self.args.write_results = self.args.output_interval > 0

        self.graph = create_countgraph(self.args)

        self.prep_results_dir()
        part_info = os.path.join(self.args.output_dir, 'partitioner-details.csv')
        self.partitioner = ConditionalPartitioner(self.graph, 
                                                  tag_density=self.args.tag_density,
                                                  info_filename=part_info)

    def parse_args(self, args):
        parser = build_counting_args(descr='Partition a sample',
                                     citations=['counting', 'SeqAn'])
        parser.add_argument('--output-dir', default='partitioned')
        parser.add_argument('samples', nargs='+')
        parser.add_argument('--save', action='store_true', default=False)
        parser.add_argument('--pairing-mode', 
                            choices=['split', 'interleaved', 'single'],
                            default='split')
        parser.add_argument('-Z', dest='norm', default=10, type=int)
        parser.add_argument('--output-interval', default=0, type=int)
        parser.add_argument('--tag-density', default=None, type=int)
        parser.add_argument('--coverage', default=10, type=int)
        parser.add_argument('--coverage-ceiling', default=None, type=int)
        
        return parser.parse_args(args)

    def write_results(self, folder, n, new_kmers):
        filename = os.path.join(folder, '{0}.csv'.format(n))
        print('# {0}: {1} tags, {2} components.'.format(n, self.partitioner.n_tags, 
                                                        self.partitioner.n_components))
        print('  writing results to file -> {0}'.format(filename))
        self.partitioner.write_components(filename)
        with open(os.path.join(folder, 'global.csv'), 'a') as fp:
            fp.write('{0}, {1}, {2}, {3}\n'.format(n, self.partitioner.n_components,
                                                 self.partitioner.n_tags, new_kmers))

    def prep_results_dir(self):
        try:
            os.mkdir(self.args.output_dir)
        except OSError as e:
            pass

        if self.args.save:
            self.args.save = os.path.join(self.args.output_dir, 'partitioner')

    def write_meta(self, n_sequences, total_kmers, fitness_func):
        meta = {'samples': self.args.samples,
                'pairing': self.args.pairing_mode,
                'K': self.args.ksize,
                'tag-density': self.partitioner.tag_density,
                'n_sequences': n_sequences,
                'n_unique_kmers': total_kmers,
                'partitioning_function': fitness_func.__class__.__module__ + '.'\
                                         + fitness_func.__class__.__name__}
        if self.args.save:
            meta['partitioner'] = self.args.save

        with open(os.path.join(self.args.output_dir, 'meta'), 'w') as fp:
            json.dump(meta, fp, indent=4)

    def run(self):

        if self.args.pairing_mode == 'split':
            samples = list(grouper(2, self.args.samples))
            for pair in samples:
                if len(pair) != 2:
                    raise ValueError('Must have even number of samples!')
        else:
            samples = self.args.samples
        
        cdef int n
        cdef int n_sequences = 0
        cdef bool paired
        cdef Sequence first, second
        cdef int new_kmers = 0
        cdef int total_kmers = 0
        cdef int print_interval = self.args.output_interval if self.args.write_results else 10000

        cdef PartitionFunction pfunc
        if self.args.coverage_ceiling is None:
            pfunc = PartitionCoverage(coverage_cutoff=self.args.coverage,
                                      graph=self.graph, partitioner=self.partitioner)
        else:
            pfunc = PartitionCoverageSlice(floor=self.args.coverage, 
                                           ceiling=self.args.coverage_ceiling,
                                           graph=self.graph, 
                                           partitioner=self.partitioner)
 
        last = 0
        for group in samples:
            if self.args.pairing_mode == 'split':
                sample_name = '{0}.{1}'.format(group[0], group[1])
                print('== Starting ({0}) =='.format(sample_name))
                reader = SplitPairedReader(FastxParser(group[0]),
                                           FastxParser(group[1]),
                                           min_length=self.args.ksize)
            else:
                sample_name = group
                print('== Starting {0} =='.format(sample_name))
                reader = BrokenPairedReader(FastxParser(group), min_length=self.args.ksize)
            for n, paired, first, second in reader:

                if n % print_interval == 0:
                    print (n, self.partitioner.n_components, self.partitioner.n_tags)
                if self.args.write_results and n > 0 and n % self.args.output_interval == 0:
                    self.write_results(self.args.output_dir, last+n, new_kmers)
                    total_kmers += new_kmers
                    new_kmers = 0

                if paired:
                    new_kmers += self.partitioner.consume_pair(first,
                                                               second,
                                                               func=pfunc)
                else:
                    new_kmers += self.partitioner.consume(first, func=pfunc)
            last = n
            n_sequences += last
            if self.args.write_results:
                self.write_results(self.args.output_dir, last, new_kmers)
                total_kmers += new_kmers
                new_kmers = 0

        if self.args.save:
            self.partitioner.save(self.args.save)

        self.write_meta(n_sequences, total_kmers, pfunc)

        return self.partitioner



