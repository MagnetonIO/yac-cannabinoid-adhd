"""YAC Cannabinoid Compound Library for ADHD Research."""

from .compounds import YACCompoundLibrary, MolecularAnalyzer, CompoundData
from .types import (
    SMILES,
    InChI,
    InChIKey,
    DrugLikenessRule,
    BindingMode,
    SignalingPathway,
    ReceptorTarget,
    ActivationMechanism,
    CompoundStrategy,
    FunctionalGroup,
    CannabinoiFunctionalGroups,
    PhysicochemicalProperties,
    BindingAffinity,
    BiasProfile,
    ProdrugDesign,
    MultiTargetProfile,
    CompoundSpecification,
)

__version__ = "1.0.0"
__all__ = [
    # Core classes
    "YACCompoundLibrary",
    "MolecularAnalyzer",
    "CompoundData",
    # Type-safe types
    "SMILES",
    "InChI",
    "InChIKey",
    # Enums
    "DrugLikenessRule",
    "BindingMode",
    "SignalingPathway",
    "ReceptorTarget",
    "ActivationMechanism",
    "CompoundStrategy",
    # Chemical types
    "FunctionalGroup",
    "CannabinoiFunctionalGroups",
    "PhysicochemicalProperties",
    "BindingAffinity",
    "BiasProfile",
    "ProdrugDesign",
    "MultiTargetProfile",
    "CompoundSpecification",
]
