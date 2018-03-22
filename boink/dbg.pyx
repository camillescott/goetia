# boink/dbg.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

include "dbg.pyx.pxi"

def make_dBG(int K, uint64_t starting_size, int n_tables,
             str storage="BitStorage", str shifter="DefaultShifter"):
    return _make_dbg(K, starting_size, n_tables, storage, shifter)
