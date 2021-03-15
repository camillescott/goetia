import argparse

from goetia import libgoetia
from goetia.dbg import dBG
from goetia.hashing import StrandAware, FwdLemireShifter, CanLemireShifter
from goetia.parsing import iter_fastx_inputs, get_fastx_args
from goetia.storage import *
from goetia.timer import measure_time

parser = argparse.ArgumentParser()
group = get_fastx_args(parser)
group.add_argument('-i', dest='inputs', nargs='+', required=True)
args = parser.parse_args()


for storage_t in [SparseppSetStorage, PHMapStorage, BitStorage, BTreeStorage]:
    for hasher_t in [FwdLemireShifter, CanLemireShifter]:
        hasher = hasher_t(31)
        if storage_t is BitStorage:
            storage = storage_t.build(int(1e9), 4)
        else:
            storage = storage_t.build()
        graph = dBG[storage_t, hasher_t].build(storage, hasher)
        consumer = dBG[storage_t, hasher_t].Processor.build(graph, 100000)

        for sample, name in iter_fastx_inputs(args.inputs, args.pairing_mode):
            print(f'dBG type: {type(graph)}')
            with measure_time():
                consumer.process(*sample)
