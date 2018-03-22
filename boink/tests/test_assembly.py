import pytest

from boink.dbg import make_dBG
from boink.assembly import make_assembler
from boink.tests.utils import *
from boink.tests.test_dbg import dbg_type
from khmer.tests.graph_structure_fixtures import *

@pytest.fixture
def asm_type(request, ksize, dbg_type):

    def build():
        G = dbg_type()
        return make_assembler(G)

    return build


def test_assembler_type(asm_type):
    asm = asm_type()

    assert asm.storage_type == asm.Graph.storage_type
    assert asm.shifter_type == asm.Graph.shifter_type

    _, _, asm_suffix = type(asm).__name__.partition('_')
    _, _, dbg_suffix = type(asm.Graph).__name__.partition('_')

    assert asm_suffix == dbg_suffix, (asm_suffix, dbg_suffix)
