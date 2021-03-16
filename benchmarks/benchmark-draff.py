import argparse
import curio

from goetia import libgoetia
from goetia.hashing import StrandAware
from goetia.parsing import iter_fastx_inputs, get_fastx_args
from goetia.processors import AsyncSequenceProcessor
from goetia.storage import *
from goetia.timer import measure_time

parser = argparse.ArgumentParser()
group = get_fastx_args(parser)
group.add_argument('-i', dest='inputs', nargs='+', required=True)
args = parser.parse_args()

for storage_t in [HLLStorage]:
    sig_t = libgoetia.sketches.UnikmerSketch[storage_t, StrandAware]
    sig = sig_t.Sketch.build(31, 10)
    proc = AsyncSequenceProcessor (
               sig_t.Processor.build(sig, 100000),
               iter_fastx_inputs(args.inputs, args.pairing_mode, names=args.names)
           )

    with measure_time():
        print(f'Storage: {storage_t}')
        curio.run(proc.start)
   
    print(sig.to_numpy())
