"""One explicit execution-size budget; not a chemistry or suitability guarantee."""

from island.exceptions import AmberToolsInputError

DEFAULT_MAX_ATOMS = 100
PROVIDED_MAX_ATOMS = 1000
AM1BCC_MAX_ATOMS = 100


def validate_size_policy(
    max_atoms, charge_method, actual_count=None, *, error_type=AmberToolsInputError
):
    ceiling = PROVIDED_MAX_ATOMS if charge_method == "provided" else AM1BCC_MAX_ATOMS
    if (
        charge_method not in ("provided", "am1bcc")
        or type(max_atoms) is not int
        or not 1 <= max_atoms <= ceiling
    ):
        raise error_type(
            f"max_atoms must be an integer in 1..{ceiling} for charge mode {charge_method}; configured={max_atoms!r}"
        )
    if actual_count is not None and (
        type(actual_count) is not int or not 1 <= actual_count <= max_atoms
    ):
        raise error_type(
            f"AmberTools actual count={actual_count}, configured limit={max_atoms} sites, charge mode={charge_method}"
        )


def policy_record(max_atoms, charge_method, actual_count):
    validate_size_policy(max_atoms, charge_method, actual_count)
    return {
        "schema": "island_ambertools_size_policy_v1",
        "max_atoms": max_atoms,
        "actual_atoms": actual_count,
        "charge_method": charge_method,
    }
