include "dbg.pyx.pxi"

def make_dBG(int K, uint64_t starting_size, int n_tables,
             str storage="BitStorage", str shifter="DefaultShifter"):
    return _make_dbg(K, starting_size, n_tables, storage, shifter)
