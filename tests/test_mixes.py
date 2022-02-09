from typing import cast
import pytest
import re
import numpy as np

from alhambra.mixes import (
    Q_,
    Component,
    FixedConcentration,
    FixedVolume,
    Mix,
    MultiFixedConcentration,
    MultiFixedVolume,
    Strand,
    VolumeError,
    WellPos,
    Reference,
    load_reference,
    nM,
    uM,
    ureg,
)
import itertools
import pandas as pd


def test_wellpos_movement():
    "Ensure WellPos movements are correct, and fail when appropriate."

    assert WellPos("A5").next_byrow() == WellPos("A6")
    assert WellPos("A12").next_byrow() == WellPos("B1")
    assert WellPos("A12").next_bycol() == WellPos("B12")

    with pytest.raises(
        ValueError, match=r"Row I \(9\) out of bounds for plate size 96"
    ):
        WellPos("H12").next_byrow()

    with pytest.raises(ValueError, match="Column 13 out of bounds for plate size 96"):
        WellPos("H12").next_bycol()

    assert WellPos("A12", platesize=384).next_byrow() == "A13"

    assert WellPos("H14", platesize=384).next_bycol() == "I14"

    assert WellPos("A12", platesize=384) == WellPos("A12", platesize=96)

    assert WellPos("D6") == WellPos(4, 6)

    assert WellPos("D6") == "D6"

    assert str(WellPos("D6")) == "D6"

    assert repr(WellPos("D8")) == 'WellPos("D8")'

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

    allbyrow96 = [f"{r}{c}" for c in range(1, 13) for r in "ABCDEFGH"]
    for x, y in itertools.pairwise(allbyrow96):
        assert WellPos(x).next_bycol() == y


def test_component():

    assert Component("test1") != Component("test2")

    assert Component("test3") == Component("test3")

    assert Component("A", 1 * uM) == Component("A", 1000 * nM)

    assert Component("A", 1 * uM) != Component("A", 1002 * nM)

    assert Component("A") != Strand("A")

    assert Component("A") != 5


def test_component_allcomps():
    ac = Component("A", 1 * uM).all_components()

    assert len(ac) == 1
    assert ac.loc["A", "component"] == Component("A", 1 * uM)
    assert np.allclose(ac.loc["A", "concentration_nM"], 1000.0)


@pytest.fixture
def reference():
    return Reference.from_csv("tests/test_reference.csv")

def test_reference_saveload(tmp_path_factory: pytest.TempPathFactory, reference: Reference):
    sf = tmp_path_factory.mktemp("exp") / "test.csv"

    r = load_reference("tests/test_reference.csv")

    assert r == reference

    r.to_csv(sf)

    r2 = Reference.from_csv(sf)

    assert r == r2

    assert r == r2.df

def test_component_with_reference(reference: Reference):
    c = Component("comp1")
    d = c.with_reference(reference)

    assert c != d
    assert np.allclose(d.concentration, ureg.Quantity(1000.0, "nM"))

    with pytest.raises(ValueError):
        Component("comp1", ureg.Quantity(150.0, "nM")).with_reference(reference)


def test_strand_with_reference(reference: Reference):
    c = Strand("strand1")
    d = c.with_reference(reference)

    assert c != d
    assert np.allclose(d.concentration, ureg.Quantity(1000.0, "nM"))
    assert d.sequence == "AGAACC"

    with pytest.raises(ValueError):
        Strand("strand1", ureg.Quantity(150.0, "nM")).with_reference(reference)

    with pytest.raises(ValueError):
        Strand("strand1", sequence="AGCTG").with_reference(reference)


def test_with_reference_get_first(
    reference: Reference, caplog: pytest.LogCaptureFixture
):
    s = Strand("strand3").with_reference(reference)

    r1 = caplog.records[0]
    assert re.match(r"Strand %s has more than one location", r1.msg)
    assert r1.args == (
        s.name,
        [("P 2", "D7"), ("P 2", "D5"), ("P 3", "D7"), ("P 4", "D7")],
    )

    assert s == Strand("strand3", Q_(1000, nM), "GGTG", plate="P 2", well="D7")


def test_with_reference_constraints_match_plate(reference: Reference, caplog):
    s = Strand("strand3", plate="P 3").with_reference(reference)
    assert s == Strand("strand3", Q_(2000, nM), "GGTG", plate="P 3", well="D7")

    c = Component("strand3", plate="P 3").with_reference(reference)
    assert c == Component("strand3", Q_(2000, nM), plate="P 3", well="D7")


def test_with_reference_constraints_match_well(reference: Reference, caplog):
    s = Strand("strand3", well="D5").with_reference(reference)
    assert s == Strand("strand3", Q_(1000, nM), "GGTG", plate="P 2", well="D5")

    c = Component("strand3", well="D5").with_reference(reference)
    assert c == Component("strand3", Q_(1000, nM), plate="P 2", well="D5")


def test_with_reference_constraints_match_seq(reference: Reference, caplog):
    s = Strand("strand3", sequence="GGTGAGG").with_reference(reference)
    assert s == Strand("strand3", Q_(2000, nM), "GGTG AGG", plate="P 4", well="D7")


def test_a_mix(reference: Reference):
    c1 = Component("comp1")
    s1 = Strand("strand1")
    s2 = Strand("strand2")
    s3 = Strand("strand3", ureg("1000 nM"), "GGTG")

    m = Mix(
        [
            MultiFixedVolume([s1, s2, s3], ureg("10 uL"), compact_display=True),
            FixedConcentration(c1, "100 nM"),
            FixedVolume(s3, ureg("10 uL")),
        ],
        name="test",
        test_tube_name="tm1",
        fixed_total_volume=ureg("50 uL"),
        fixed_concentration="strand3",
    ).with_reference(reference)

    assert m.buffer_volume == ureg("5 uL")
    assert np.allclose(m.concentration, ureg("400 nM"))

    mdt = m._repr_markdown_().splitlines()

    assert (
        re.match(
            r"Table: Mix: test, Conc: 400.00 nM, Total Vol: 50.00 µl, Component Count: 4, Test tube name: tm1",
            mdt[0],
        )
        is not None
    )

    # FIXME: more tests of output

    # FIXME: test of chained mix


def test_multifixedconc_min_volume(reference: Reference):
    s1 = Strand("strand1", "400 nM")
    s2 = Strand("strand2", "200 nM")

    m = Mix(
        [MultiFixedConcentration([s1, s2], "50 nM", min_volume="20 uL")],
        name="test",
        fixed_total_volume="100 uL",
    )
    
    with pytest.raises(VolumeError):
        m.table()

    m.fixed_total_volume = "200 uL"  # type: ignore  # Mypy doesn't understand on_setattr

    m.table()

def test_non_plates():
    s1 = Strand("s1", "200 nM", plate="tube")

    s2 = Strand("s2", "200 nM", plate="tube")

    s3 = Strand("s3", "400 nM", plate="tube")

    s4 = Strand("s4", "400 nM", plate="a different tube")

    m = Mix([MultiFixedVolume([s1, s2, s3, s4], "1 uL", equal_conc="min_volume")], "test")

    m.table()

    ml = m.mixlines()

    assert len(ml) == 4