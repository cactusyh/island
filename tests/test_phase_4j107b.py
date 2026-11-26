"""J107B source-bound chemistry coverage and GAFF graph checks."""

import pytest


def test_profile_v2_installs_without_mutating_base(tmp_path):
    from island.forcefields.emc_profiles import install

    destination = tmp_path / "field"
    manifest = install("opls-polymer-v2", destination)
    assert manifest["numerical_parameters_changed"] is False
    assert set(manifest["source_residues"]) == {
        "ETCL",
        "ETBR",
        "ETIO",
        "MEAN",
        "DMAN",
        "TMAN",
        "ANIS",
    }
    assert (destination / "profile.json").exists()


def test_source_and_heldout_profile_checks(tmp_path):
    from island.forcefields.emc_profiles import install
    from scripts.validate_j107b_coverage import check

    field = tmp_path / "field"
    install("opls-polymer-v2", field)
    result = check(field, tmp_path / "checks")
    assert result["passed"]
    assert result["source_models"] == 8 and result["development_models"] == 5


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
def test_common_gaff_graph_coverage(family, tmp_path):
    from scripts.validate_gaff_graph_coverage import check

    result = check(family)
    assert result["attempted"] == 12
    assert result["assigned"] == result["attempted"]
    assert all(r["native_match"] for r in result["records"])


def test_profile_unknown_name_rejected(tmp_path):
    from island.forcefields.emc_profiles import install

    with pytest.raises(ValueError, match="Unknown EMC profile"):
        install("opls-polymer-unknown", tmp_path / "field")
