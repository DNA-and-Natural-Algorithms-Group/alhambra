import pytest
from alhambra import *

try:
    import xgrow

    xgrow_available = True
except ImportError:
    xgrow_available = False


def test_xor_example():
    ts = TileSet.from_file("./examples/xor_ribbon.yaml")

    assert ts.tiles["t01_1"].edges["N"].name == "r0_v1*"


# @pytest.mark.skipif(not xgrow_available, reason="xgrow not available")
# def test_xor_example_xgrow():
#     ts = TileSet.from_file("./examples/xor_ribbon.yaml")

#     out = ts.run_xgrow(smax=80, window=False)

#     ts.create_abstract_diagram(out)
