"""Install explicitly selected typing profiles on a pinned EMC distribution."""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory

from .emc import _field_sources

PCFF = {
    "cff/pcff/pcff.frc": "22974f591c98016c62e6986d8ca831713dbde38a203ed0e5805409b17fb00013",
    "cff/pcff/pcff_templates.dat": "c43b923bab8e68cf534902c79e54de01b502e6e25dd653910ba2aa41f06c45e6",
}


def carboxamide(template):
    # Both the carbonyl carbon and the amide N have explicit native PCFF types.
    rule = """type: n
    ! neutral carboxamide nitrogen (carbonyl substituent is C)
    template: (>N(-C(=O)(-C))(-*)(-*))
    atom_test:1
      formal_charge:0
      aromaticity:NON_AROMATIC
    end_test
end_type

"""
    return template.replace("precedence:", rule + "precedence:").replace(
        "(na(n3m)(n4m)(n_2))", "(na(n)(n3m)(n4m)(n_2))"
    )


OPLS = {
    "opls/2024/opls-aa.define": "ad94996fb33f437053b46c42af192987e435fdfdad14106a37c216805e48a922",
    "opls/2024/opls-aa.prm": "a0becae1a8f4d9faf9a419d8cc626ce44c2fb6b364838607dc1aa11c8aeba469",
    "opls/2024/opls-aa.top": "f46d83cc2da20d803a9f24f0aa5ac2ca1393d027d20e845c352d7516c42acbc1",
}


def epoxy_amine(template):
    # Native types from EOXI, ANIS, EMET and the primary/secondary/tertiary
    # aliphatic amine residues. Extend graph matching to substituted rings.
    rules = """1000	c4oe	C	-	-	0	^3C(^3C)(^3O)(C)(H)
1001	c4oe	C	-	-	0	C(O(^6c))(C)(H)(H)
1002	c4oe	C	-	-	0	C(O(^6c))(^3C)(H)(H)
1003	c4n	C	-	-	0	^6C(^6C)(^6C)(N)(H)
1004	c4n	C	-	-	0	C(N)(C)(H)(H)
1005	h1n	H	-	-	0	H(N(^6C)(H))
1006	h1n	H	-	-	0	H(N(^6C)(C))
1007	h1n	H	-	-	0	H(N(C)(H))
1008	h1cn	H	-	-	0	H(^6C(^6C)(^6C)(N))
1009	h1cn	H	-	-	0	H(C(N)(C)(H))
1010	c4oe	C	-	-	0	C(O(C))(C)(H)(H)
"""
    end = "ITEM\tEND\t# RULES"
    return template.replace(end, rules + end)


def polymer_opls(template):
    """Extend native ester/ether and N-alkyl amide contexts from source residues."""
    template = epoxy_amine(template)
    # MACE's methyl ester H rule must retain its carbonyl environment.
    # The original short pattern also matches the different EMET/DMET H type.
    template = template.replace(
        "256\th1co\tH\tMACE\tH31\t0\tH(C(O)(H)(H))",
        "256\th1co\tH\tMACE\tH31\t0\tH(C(O(C(=O)(*)))(H)(H))",
    )
    patterns = [
        ("c4oe", "C", "C(O(C))(C)(C)(H)"),
        ("c4oe", "C", "C(O(C))(C)(C)(C)"),
        ("c4oe", "C", "C(O(C))(H)(H)(H)"),
        ("h1co", "H", "H(C(C)(C)(O(C(=O)(*))))"),
        ("c4n", "C", "C(N(C(C)(H)(H))(H))(H)(H)(H)"),
    ]
    for other in ("H", "C"):
        nitrogen = f"N(C(=O)(*))({other})"
        for substituents in ("(H)(H)(H)", "(C)(H)(H)", "(C)(C)(H)", "(C)(C)(C)"):
            patterns.append(("c4", "C", f"C({nitrogen}){substituents}"))
        for substituents in ("(H)(H)", "(C)(H)", "(C)(C)"):
            patterns.append(("h1", "H", f"H(C({nitrogen}){substituents})"))
    rules = "".join(
        f"{1100 + i}\t{kind}\t{element}\t-\t-\t0\t{pattern}\n"
        for i, (kind, element, pattern) in enumerate(patterns)
    )
    return template.replace("ITEM\tEND\t# RULES", rules + "ITEM\tEND\t# RULES")


def polymer_opls_v2(template):
    """Source-defined alkyl halide, neutral amine and aryl ether contexts.

    Numerical parameters remain byte-identical. ETCL/ETBR/ETIO, MEAN/DMAN/
    TMAN and ANIS supply the native types and increments; their MASS entries
    explicitly identify the corresponding chemical classes.
    """
    template = polymer_opls(template)
    # Two stock generic c4 rules override chloro-alkane typing at substituted
    # carbons but c4/cl1 has no increment. Use the declared chloro-alkane type.
    for number in (91, 97):
        old = f"{number}\tc4\tC\t-\t-\t0\t"
        if template.count(old) != 1:
            raise ValueError("Chloro-alkane source rule changed")
        template = template.replace(old, f"{number}\tc4cl\tC\t-\t-\t0\t")
    patterns = []
    for atom, kind in (("[Br]", "c4br"), ("I", "c4i")):
        for tails in ("(C)(C)(H)", "(C)(C)(C)"):
            patterns.append((kind, "C", f"C({atom}){tails}"))
    # All carbon branches attached to N are saturated, with four neighbors;
    # this excludes carbonyl C, aromatic c, nitriles and charged N contexts.
    for nitrogen in (
        "N(H)(H)",
        "N(C(~*)(~*)(~*))(H)",
        "N(C(~*)(~*)(~*))(C(~*)(~*)(~*))",
    ):
        for tails in ("(H)(H)(H)", "(C)(H)(H)", "(C)(C)(H)", "(C)(C)(C)"):
            patterns.append(("c4n", "C", f"C({nitrogen}){tails}"))
        for tails in ("(H)(H)", "(C)(H)", "(C)(C)"):
            patterns.append(("h1cn", "H", f"H(C({nitrogen}){tails})"))
    for tails in ("(C)(C)(H)", "(C)(C)(C)"):
        patterns.append(("c4oe", "C", f"C(O(^6c)){tails}"))
    rules = "".join(
        f"{1200 + i}\t{kind}\t{element}\t-\t-\t0\t{pattern}\n"
        for i, (kind, element, pattern) in enumerate(patterns)
    )
    return template.replace("ITEM\tEND\t# RULES", rules + "ITEM\tEND\t# RULES")


def install(profile, directory, *, emc_root=None):
    if profile == "opls-polymer-v2":
        field, pins, edits = (
            "opls/2024/opls-aa",
            OPLS,
            {"opls/2024/opls-aa.top": polymer_opls_v2},
        )
    elif profile == "opls-polymer-v1":
        field, pins, edits = (
            "opls/2024/opls-aa",
            OPLS,
            {"opls/2024/opls-aa.top": polymer_opls},
        )
    elif profile == "opls-epoxy-amine-v1":
        field, pins, edits = (
            "opls/2024/opls-aa",
            OPLS,
            {"opls/2024/opls-aa.top": epoxy_amine},
        )
    elif profile == "pcff-carboxamide-v1":
        field, pins, edits = "pcff", PCFF, {"cff/pcff/pcff_templates.dat": carboxamide}
    else:
        raise ValueError(f"Unknown EMC profile: {profile}")
    _, _, source, hashes = _field_sources(field, emc_root, None)
    if hashes != pins:
        raise ValueError("Installed field differs from the validated source")
    out = Path(directory)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise FileExistsError(out)
    with TemporaryDirectory(dir=out.parent) as tmp:
        stage = Path(tmp) / "fields"
        stage.mkdir()
        result = {}
        for name in hashes:
            path = stage / name
            path.parent.mkdir(parents=True, exist_ok=True)
            raw = (
                edits[name]((source / name).read_text()).encode()
                if name in edits
                else (source / name).read_bytes()
            )
            path.write_bytes(raw)
            result[name] = sha256(raw).hexdigest()
        manifest = {"profile": profile, "base": hashes, "installed": result}
        if profile == "opls-polymer-v2":
            manifest["source_residues"] = [
                "ETCL",
                "ETBR",
                "ETIO",
                "MEAN",
                "DMAN",
                "TMAN",
                "ANIS",
            ]
            manifest["numerical_parameters_changed"] = False
        (stage / "profile.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(out)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile")
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(install(args.profile, args.directory), indent=2))
