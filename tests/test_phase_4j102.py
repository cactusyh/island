from dataclasses import asdict

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from island.epoxy.chemistry import AMINE, THIOL, react, reference_reactions
from island.epoxy.cure import CureConfig
from island.epoxy.network import select


def test_thiol_opens_epoxide_once_and_transfers_the_same_h():
    _, hardener, [(pre, event, post)] = reference_reactions(template=THIOL)
    assert len(THIOL.donors(hardener)) == 4
    assert set(asdict(event)) == {"carbon", "oxygen", "sulfur", "hydrogen"}
    assert THIOL.event(asdict(event)) == event
    assert pre.GetBondBetweenAtoms(event.sulfur, event.hydrogen)
    assert post.GetBondBetweenAtoms(event.oxygen, event.hydrogen)
    assert post.GetBondBetweenAtoms(event.carbon, event.sulfur)
    assert not post.GetBondBetweenAtoms(event.carbon, event.oxygen)
    assert not THIOL.donors(post)[event.sulfur]
    assert CalcMolFormula(pre) == CalcMolFormula(post)
    assert Chem.GetFormalCharge(pre) == Chem.GetFormalCharge(post) == 0
    for mol, events, template in [
        (pre, [event, event], THIOL),
        (post, [event], THIOL),
        (pre, [event], AMINE),
    ]:
        with pytest.raises(ValueError):
            react(mol, events, template)


def test_chemistry_selection_is_unambiguous():
    assert CureConfig().pacm == 4 and CureConfig().template == AMINE
    assert CureConfig(petmp=2).template == THIOL
    with pytest.raises(ValueError):
        CureConfig(petmp=2, pacm=2)
    with pytest.raises(ValueError):
        CureConfig(pacm=0)


def test_periodic_sulfur_angle_rejects_linear_approach():
    _, _, [(mol, event, _)] = reference_reactions(template=THIOL)
    xyz = np.full((mol.GetNumAtoms(), 3), 10.0)
    sulfur = event.sulfur
    neighbor = next(
        a.GetIdx()
        for a in mol.GetAtomWithIdx(sulfur).GetNeighbors()
        if a.GetAtomicNum() == 6
    )
    xyz[list(THIOL.donors(mol))] = [10, 15, 10]
    xyz[sulfur], xyz[event.carbon] = [19, 4, 5], [1, 4, 5]
    box = np.full(3, 20.0)

    def choose(pos):
        return select(
            mol,
            pos,
            box,
            3,
            np.random.default_rng(3),
            1,
            template=THIOL,
            donor_angle_tolerance=45,
        )

    xyz[neighbor] = [18, 4, 5]
    events, record = choose(xyz)
    assert not events and record["rejected"]["donor_angle"] == 1
    xyz[neighbor] = [19, 5, 5]
    for pos in (xyz, (xyz + [3, 6, 9]) % box):
        events, record = choose(pos)
        assert events == [event]
        assert record["events"][0]["donor_angle_deg"] == pytest.approx(90)
