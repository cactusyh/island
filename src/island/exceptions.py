"""Package-specific exceptions."""


class IslandError(Exception):
    """Base exception for ISLAND."""


class TopologyError(IslandError):
    """Raised for invalid chemical topology operations."""


class CoordinateError(IslandError):
    """Raised for invalid coordinate operations."""


class ValidationError(IslandError):
    """Raised when a model fails consistency validation."""


class ReactionError(IslandError):
    """Raised for invalid reaction transformation operations."""


class ForceFieldError(IslandError):
    """Raised for force-field parameterization failures."""


class AtomTypingError(ForceFieldError):
    """Base exception for atom-typing operations."""


class InvalidAtomTypingRuleError(AtomTypingError):
    """Raised when an atom-typing rule or ruleset is invalid."""


class UnsupportedAtomTypingError(AtomTypingError):
    """Raised when an atom-typing engine cannot type the supplied chemistry."""


class IncompleteAtomTypingError(AtomTypingError):
    """Raised when strict atom typing leaves untyped or ambiguous sites."""

    def __init__(self, message: str, *, result: object) -> None:
        super().__init__(message)
        self.result = result


class ParameterAssignmentError(ForceFieldError):
    """Base exception for numerical parameter-assignment operations."""


class InvalidParameterDefinitionError(ParameterAssignmentError):
    """Raised when a parameter record or library is invalid."""


class InvalidTypingResultError(ParameterAssignmentError):
    """Raised when atom-typing input is stale or structurally inconsistent."""


class InvalidParameterAssignmentResultError(ParameterAssignmentError):
    """Raised when a parameter-assignment result is internally inconsistent."""


class UnsupportedParameterRequirementError(ParameterAssignmentError):
    """Raised when required interactions are outside the supported scope."""


class IncompleteParameterAssignmentError(ParameterAssignmentError):
    """Raised when strict parameter assignment is incomplete."""

    def __init__(self, message: str, *, result: object) -> None:
        super().__init__(message)
        self.result = result


class ChargeAssignmentError(ForceFieldError):
    """Base exception for partial-charge assignment operations."""


class InvalidChargeDefinitionError(ChargeAssignmentError):
    """Raised for malformed provided charges or charge-table definitions."""


class InvalidChargeAssignmentResultError(ChargeAssignmentError):
    """Raised when a charge-assignment result is internally inconsistent."""


class IncompleteChargeAssignmentError(ChargeAssignmentError):
    """Raised when strict charge assignment is incomplete or charge-inconsistent."""

    def __init__(self, message: str, *, result: object) -> None:
        super().__init__(message)
        self.result = result


class NonbondedPolicyError(ForceFieldError):
    """Raised for unsupported or invalid nonbonded-policy operations."""


class ParameterCompositionError(ForceFieldError):
    """Raised when typed parameters, charges, and policy cannot be composed."""


class AmberImportError(ForceFieldError):
    """Base error for importing an already resolved Amber topology."""


class UnsupportedAmberFeatureError(AmberImportError):
    """The source topology cannot be represented faithfully by this adapter."""


class InvalidAmberMappingError(AmberImportError):
    """Source atom indices do not map bijectively to authoritative sites."""


class InvalidAmberImportResultError(AmberImportError):
    """An imported result is stale or internally inconsistent."""


class ChemistryError(IslandError):
    """Base exception for optional chemistry-layer operations."""


class RDKitConversionError(ChemistryError):
    """Raised when conversion to or from an RDKit molecule fails."""


class UnsupportedRepresentationError(RDKitConversionError):
    """Raised when a representation cannot be expressed by an adapter."""


class MissingConformerError(RDKitConversionError):
    """Raised when an RDKit molecule has no required conformer."""


class EmbeddingError(ChemistryError):
    """Raised when three-dimensional coordinate embedding fails."""


class PSMILESError(ChemistryError):
    """Base exception for PSMILES parsing and interpretation."""


class InvalidRepeatUnitError(PSMILESError):
    """Raised when a PSMILES repeat unit is invalid or unsupported."""


class PolymerBuildError(ChemistryError):
    """Raised when a polymer molecular graph cannot be constructed."""


class StereochemistryError(ChemistryError):
    """Base exception for chemical stereochemistry operations."""


class TacticityError(StereochemistryError):
    """Raised for invalid or inapplicable tacticity requests."""


class UnsupportedStereochemistryError(TacticityError):
    """Raised when requested stereochemistry exceeds the supported model."""


class ConformationError(IslandError):
    """Base exception for coordinate-only conformation operations."""


class UnsupportedConformationError(ConformationError):
    """Raised when a conformation generator cannot handle a molecular system."""


class ConformationGenerationError(ConformationError):
    """Raised when conformation generation exhausts its retry policy."""

    def __init__(
        self,
        message: str,
        *,
        repeat_index: int | None = None,
        attempts: int = 0,
        rejected_trials: int = 0,
        rollback_count: int = 0,
    ) -> None:
        super().__init__(message)
        self.repeat_index = repeat_index
        self.attempts = attempts
        self.rejected_trials = rejected_trials
        self.rollback_count = rollback_count
