from black import re
import numpy as np
import pytest

from alhambra.mixes import MOLARITY_TOLERANCE, Component, FixedConcentration, FixedVolume, Mix, MultiFixedVolume, Strand, WellPos, load_reference, nM, uM, UR
import itertools
import pandas as pd

def test_wellpos_movement():
    "Ensure WellPos movements are correct, and fail when appropriate."

    assert WellPos("A5").next_byrow() == WellPos("A6")
    assert WellPos("A12").next_byrow() == WellPos("B1")
    assert WellPos("A12").next_bycol() == WellPos("B12")

    with pytest.raises(ValueError, match=r"Row I \(9\) out of bounds for plate size 96"):
        WellPos("H12").next_byrow()

    with pytest.raises(ValueError, match="Column 13 out of bounds for plate size 96"):
        WellPos("H12").next_bycol()

    assert WellPos("A12", platesize=384).next_byrow() == "A13"

    assert WellPos("H14", platesize = 384).next_bycol() == "I14"

    assert WellPos("A12", platesize = 384) == WellPos("A12", platesize = 96)

    assert WellPos("D6") == WellPos(4, 6)

    assert WellPos("D6") == "D6"

    assert str(WellPos("D6")) == "D6"


    assert repr(WellPos("D8")) == "WellPos(\"D8\")"

    assert WellPos("C8").key_byrow() == (3, 8)

    assert WellPos("C8").key_bycol() == (8, 3)

    assert WellPos("D8") == WellPos(WellPos("D8"))

def test_invalid_wellrefs():
    with pytest.raises(ValueError):
        WellPos("A14")

    with pytest.raises(ValueError):
        WellPos("Q14", platesize=384)

    with pytest.raises(ValueError):
        WellPos("H25", platesize=384)

    with pytest.raises(ValueError, match="Plate size 1536 not supported"):
        WellPos("A1", platesize=1536)

    assert WellPos("D8") != str

    with pytest.raises(TypeError):
        WellPos(5.3)

    with pytest.raises(ValueError):
        WellPos("i123nvalid string")


def test_all_wellref_96():
    allbyrow96 = [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)]
    for x, y in itertools.pairwise(allbyrow96):
        assert WellPos(x).next_byrow() == y


    allbyrow96 = [f"{r}{c}"  for c in range(1, 13) for r in "ABCDEFGH"]
    for x, y in itertools.pairwise(allbyrow96):
        assert WellPos(x).next_bycol() == y

def test_component():

    assert Component("test1") != Component("test2")

    assert Component("test3") == Component("test3")

    assert Component("A", 1*uM) == Component("A", 1000*nM)

    assert Component("A", 1*uM) != Component("A", 1002*nM)

    assert Component("A") != Strand("A")

    assert Component("A") != 5

def test_component_allcomps():
    ac = Component("A", 1*uM).all_components()

    assert len(ac) == 1
    assert ac.loc["A", "component"] == Component("A", 1*uM)
    assert np.allclose(ac.loc["A", "concentration_nM"], 1000.0)
@pytest.fixture
def reference():
    return load_reference("tests/test_reference.csv")

def test_component_with_reference(reference: pd.DataFrame):
    c = Component("comp1")
    d = c.with_reference(reference)

    assert c != d
    assert abs(d.concentration - UR.Quantity(1000.0, "nM")) < MOLARITY_TOLERANCE

    with pytest.raises(ValueError):
        Component("comp1", UR.Quantity(150.0, "nM")).with_reference(reference)

def test_strand_with_reference(reference: pd.DataFrame):
    c = Strand("strand1")
    d = c.with_reference(reference)

    assert c != d
    assert abs(d.concentration - UR.Quantity(1000.0, "nM")) < MOLARITY_TOLERANCE
    assert d.sequence == "AGAACC"

    with pytest.raises(ValueError):
        Strand("strand1", UR.Quantity(150.0, "nM")).with_reference(reference)

    with pytest.raises(ValueError):
        Strand("strand1", sequence="AGCTG").with_reference(reference)

def test_a_mix(reference: pd.DataFrame):
    c1 = Component("comp1")
    s1 = Strand("strand1")
    s2 = Strand("strand2")
    s3 = Strand("strand3", UR("1000 nM"), "GGTG")

    m = Mix("test",
    [MultiFixedVolume([s1, s2, s3], UR("10 uL"), compact_display=True),
     FixedConcentration(c1, UR("100 nM")),
     FixedVolume(s3, UR("10 uL"))
    ], UR("50 uL"), fixed_concentration='strand3'
    ).with_reference(reference)

    assert m.buffer_volume == UR("5 uL")
    assert m.concentration == UR("400 nM")

    mdt = m._repr_markdown_().splitlines()

    assert re.match(r"Table: Mix: test, Conc: 400.00 nM, Total Vol: 50.00", mdt[0]) is not None

    # FIXME: more tests of output

    # FIXME: test of chained mix