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
