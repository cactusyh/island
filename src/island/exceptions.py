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
