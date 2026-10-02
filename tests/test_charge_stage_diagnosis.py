"""Independent decimal sums and explicit lineage failures for diagnostic readers."""

from decimal import Decimal

import pytest

from scripts.diagnose_oligomer_residuals import ac, charge, inventory, mol2, prmtop, sqm


def test_printed_sqm_sum_is_not_printed_total_and_maps_named_input():
    expected = {"C019": "C", "H003": "H"}
    q, total = sqm(
        "Atom    Element       Mulliken Charge\n1 H 0.073\n2 C -0.072\n Total Mulliken Charge = 0.000",
        " 1 H003 0 0 0\n 6 C019 1 0 0",
        expected,
    )
    assert q == {"H003": Decimal(".073"), "C019": Decimal("-.072")}
    assert sum(q.values()) == Decimal(".001")
    assert Decimal(total) == 0
    with pytest.raises(ValueError, match="element"):
        sqm(
            "Element       Mulliken Charge\n1 C 0.073\n2 C -0.072\nTotal Mulliken Charge = 0.000",
            " 1 H003 0 0 0\n 6 C019 1 0 0",
            expected,
        )


def test_decimal_ac_mol2_and_prmtop_precision():
    expected = {"C019": "C"}
    assert ac("ATOM 1 C019MOL 1 0 0 0 -0.158500 c3", expected)["C019"] == Decimal(
        "-.1585"
    )
    assert mol2(
        "@<TRIPOS>ATOM\n1 C019 0 0 0 c3 1 MOL -0.158500\n@<TRIPOS>BOND", expected
    )["C019"] == Decimal("-.1585")
    text = "%FLAG ATOM_NAME\n%FORMAT(20a4)\nC019\n%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n       6\n%FLAG CHARGE\n%FORMAT(5E16.8)\n  1.82223000E+00\n"
    assert prmtop(text, expected)["C019"] == Decimal(".1")


@pytest.mark.parametrize("rows", [[("C019", "0"), ("C019", "0")], [("H000", "0")], []])
def test_exact_inventory(rows):
    with pytest.raises(ValueError):
        inventory(rows, {"C019": "C"})


@pytest.mark.parametrize("token", ["NaN", "Infinity", "-Infinity"])
def test_nonfinite_charge(token):
    with pytest.raises(ValueError):
        charge(token)
