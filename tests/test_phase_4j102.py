from dataclasses import asdict

import pytest

pytest.importorskip("scipy")
pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from island.epoxy.chemistry import AMINE, THIOL, react, reference_reactions
from island.epoxy.cure import CureConfig


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
