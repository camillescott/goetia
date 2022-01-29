from goetia import libgoetia


typenames = [(libgoetia.FwdLemireShifter, 'FwdLemireShifter'),
             (libgoetia.CanLemireShifter, 'CanLemireShifter'),
             (libgoetia.FwdUnikmerShifter, 'FwdUnikmerShifter'),
             (libgoetia.CanUnikmerShifter, 'CanUnikmerShifter')]
types = [_type for _type, _name in typenames]


UKHS = libgoetia.UKHS
HashExtender = libgoetia.HashExtender
extender_selector_t = libgoetia.extender_selector_t

for hasher_t, name in typenames:
    globals()[name] = hasher_t


Canonical = libgoetia.Canonical['uint64_t']
StrandAware  = libgoetia.Hash['uint64_t']
