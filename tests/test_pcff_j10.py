"""J10 source accounting and graph-domain controls; real parity is separate."""

import pytest

from island.builders import build_linear_polymer
from island.exceptions import PCFFError
from island.forcefields.pcff import select_pcff_profile
from island.forcefields.pcff.registry import FUSED_BENZENOID, LINKED_BENZENOID
from island.forcefields.pcff.source import PIN


def test_final_fused_polymer_domain_is_not_old_bounded_profile():
    s = build_linear_polymer("[*:1]c1ccc2cc([*:2])ccc2c1", dp=2, generate_3d=False)
    before = s.to_dict()
    select_pcff_profile(
        FUSED_BENZENOID, sha256=PIN["sha256"]
    ).profile().validate_system(s)
    with pytest.raises(PCFFError):
        select_pcff_profile(
            LINKED_BENZENOID, sha256=PIN["sha256"]
        ).profile().validate_system(s)
    assert s.to_dict() == before


@pytest.mark.parametrize(
    "psmiles", ["[*:1]c1ccnc([*:2])c1", "[*:1]c1cc([*:2])oc1", "[*:1]CC[*:2]"]
)
def test_fused_domain_near_misses(psmiles):
    s = build_linear_polymer(psmiles, dp=1, generate_3d=False)
    with pytest.raises(PCFFError):
        select_pcff_profile(
            FUSED_BENZENOID, sha256=PIN["sha256"]
        ).profile().validate_system(s)


def test_source_row_inventory_never_promotes_bounded_observations(
    tmp_path, monkeypatch
):
    from test_pcff_operational_profile import synthetic_profile_source

    from island.forcefields.pcff import load_pcff_source
    from island.forcefields.pcff.row_coverage import pcff_source_row_ledger

    # The source fixture is explicitly synthetic and override-isolated.
    path, _ = synthetic_profile_source.__wrapped__(tmp_path, monkeypatch)
    source = load_pcff_source(path)
    data = pcff_source_row_ledger(source)
    assert data["rows"] and not data["full_source_complete"]
    rid = data["rows"][0]["id"]
    checked = pcff_source_row_ledger(source, verified_rows={rid: ["synthetic_case"]})
    assert checked["rows"][0]["bounded_verified_cases"] == ["synthetic_case"]
    assert not checked["rows"][0]["globally_verified"]
    with pytest.raises(PCFFError):
        pcff_source_row_ledger(source, verified_rows={"nonexistent": ["case"]})
