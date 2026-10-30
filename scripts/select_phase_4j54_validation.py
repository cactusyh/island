"""Retain a raw historical AmberTools failure and its exact-node successful retry.

The full-suite union contains every collected node once. The failed original XML
and log remain unchanged; the passed subset is explicitly an excerpt, not a new
pytest run. Only the named external-tool failure may be replaced by its retry.
"""

import argparse
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

from scripts.phase_4j54 import signed, write
from scripts.validate_phase_4j47_receipt import pin, require

FAILED_NODE = "tests/test_phase_4j31_execution.py::test_real_final_graph_source_completeness_and_exact_execution[crosslinked-gaff2]"


def node(case):
    return case.attrib["classname"].replace(".", "/") + ".py::" + case.attrib["name"]


def selection(root, *, create=False):
    root = Path(root)
    raw = root / "validation/full-5.xml"
    retry = root / "validation/full-5-rerun.xml"
    tree = ET.parse(raw)
    cases = tree.findall(".//testcase")
    failed = [
        c for c in cases if c.find("failure") is not None or c.find("error") is not None
    ]
    require(
        len(failed) == 1 and node(failed[0]) == FAILED_NODE,
        "Unexpected initial-suite failure; cannot substitute retry",
    )
    require(
        "Charge assignment is incomplete or inconsistent with component formal charges"
        in ET.tostring(failed[0], encoding="unicode"),
        "Changed historical external-tool diagnostic",
    )
    retried = ET.parse(retry).findall(".//testcase")
    require(
        len(retried) == 1
        and node(retried[0]) == FAILED_NODE
        and retried[0].find("failure") is None
        and retried[0].find("error") is None
        and retried[0].find("skipped") is None,
        "Retry does not independently pass the exact node",
    )
    kept = [c for c in cases if c not in failed]
    excerpt = root / "validation/full-5-passed.xml"
    if create:
        outer = ET.Element("testsuites")
        suite = ET.SubElement(
            outer,
            "testsuite",
            {
                "name": "passed excerpt of raw full-5.xml; not a new pytest execution",
                "tests": str(len(kept)),
                "failures": "0",
                "errors": "0",
                "skipped": str(sum(c.find("skipped") is not None for c in kept)),
            },
        )
        suite.extend(deepcopy(kept))
        ET.ElementTree(outer).write(excerpt, encoding="utf-8", xml_declaration=True)
    selected = ET.parse(excerpt).findall(".//testcase")
    require(
        [ET.tostring(c) for c in kept] == [ET.tostring(c) for c in selected],
        "Passed excerpt changed original successful nodes",
    )
    result = signed(
        {
            "initial_xml_sha256": pin(raw),
            "retry_xml_sha256": pin(retry),
            "passed_excerpt_sha256": pin(excerpt),
            "initial_failure_node": FAILED_NODE,
            "initial_failure_diagnostic": "AmberTools charge assignment incomplete/inconsistent in historical J31 GAFF2 test",
            "retry_status": "passed",
            "historical_code_changed": False,
            "charge_or_tolerance_change": False,
            "selected_nodes": len(selected) + 1,
            "selection": "exact-node retry union; original failed run retained",
        }
    )
    path = root / "validation/retry-selection.json"
    if create:
        write(path, result)
    else:
        require(
            json.loads(path.read_text()) == result, "Changed retry selection record"
        )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--create", action="store_true")
    a = p.parse_args()
    print(json.dumps(selection(a.directory, create=a.create), sort_keys=True))


if __name__ == "__main__":
    main()
