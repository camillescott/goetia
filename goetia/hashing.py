from goetia import libgoetia


typenames = [(libgoetia.hashing.FwdLemireShifter, 'FwdLemireShifter'),
             (libgoetia.hashing.CanLemireShifter, 'CanLemireShifter'),
             (libgoetia.hashing.FwdUnikmerShifter, 'FwdUnikmerShifter'),
             (libgoetia.hashing.CanUnikmerShifter, 'CanUnikmerShifter')]
types = [_type for _type, _name in typenames]


UKHS = libgoetia.hashing.UKHS
HashExtender = libgoetia.hashing.HashExtender
extender_selector_t = libgoetia.hashing.extender_selector_t

for hasher_t, name in typenames:
    globals()[name] = hasher_t


Canonical = libgoetia.hashing.Canonical['uint64_t']
StrandAware  = libgoetia.hashing.Hash['uint64_t']