#!/usr/bin/env python3
"""
Type-safe chemical type definitions for YAC Compound Library.

Provides strongly-typed representations of molecular structures,
functional groups, and pharmacological properties.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Literal, NewType
from abc import ABC, abstractmethod


# Type-safe SMILES representation
SMILES = NewType("SMILES", str)
InChI = NewType("InChI", str)
InChIKey = NewType("InChIKey", str)


class DrugLikenessRule(Enum):
    """Drug-likeness evaluation rules."""
    LIPINSKI = auto()  # Rule of Five
    VEBER = auto()  # Oral bioavailability
    GHOSE = auto()  # Lead-likeness
    EGAN = auto()  # BBB permeability


class BindingMode(Enum):
    """Receptor binding modes."""
    ORTHOSTERIC = auto()  # Primary binding site
    ALLOSTERIC = auto()  # Secondary modulator site
    BITOPIC = auto()  # Spans both sites


class SignalingPathway(Enum):
    """G protein-coupled receptor signaling pathways."""
    Gi = auto()  # Inhibitory G protein
    Gs = auto()  # Stimulatory G protein
    Gq = auto()  # Phospholipase C pathway
    BETA_ARRESTIN = auto()  # Arrestin-mediated signaling
    ERK = auto()  # Extracellular signal-regulated kinase


class ReceptorTarget(Enum):
    """Pharmacological receptor targets."""
    CB1 = "Cannabinoid Receptor 1"
    CB2 = "Cannabinoid Receptor 2"
    GHSR1A = "Ghrelin Receptor"
    MC4R = "Melanocortin-4 Receptor"
    NPY_Y1 = "Neuropeptide Y Receptor Y1"
    OX1R = "Orexin Receptor 1"
    OX2R = "Orexin Receptor 2"


class ActivationMechanism(Enum):
    """Prodrug activation mechanisms."""
    CYP2D6 = "Cytochrome P450 2D6"
    CYP3A4 = "Cytochrome P450 3A4"
    CES1 = "Carboxylesterase 1"
    CES2 = "Carboxylesterase 2"
    ESTERASE = "Non-specific esterase"
    ACID_LABILE = "pH-dependent hydrolysis"


class CompoundStrategy(Enum):
    """Drug design strategy categories."""
    BIASED_SIGNALING = "A"
    PRODRUG_ACTIVATION = "B"
    MULTI_TARGET = "C"


@dataclass(frozen=True)
class FunctionalGroup:
    """Immutable representation of a chemical functional group."""
    name: str
    smarts: str
    description: str

    def __post_init__(self):
        if not self.smarts:
            raise ValueError("SMARTS pattern cannot be empty")


# Common cannabinoid functional groups
class CannabinoiFunctionalGroups:
    """Type-safe functional group definitions for cannabinoid scaffolds."""

    PHENOLIC_OH = FunctionalGroup(
        name="C1-OH",
        smarts="[OH]c1ccccc1",
        description="Phenolic hydroxyl at C1 position (required for CB1 binding)"
    )

    METHYL_ETHER = FunctionalGroup(
        name="C1-OCH3",
        smarts="[CH3]Oc1ccccc1",
        description="Methyl ether at C1 (reduces efficacy to ~50%)"
    )

    DIFLUOROMETHYL = FunctionalGroup(
        name="C9-CHF2",
        smarts="[CH](F)F",
        description="Difluoromethyl at C9 (alters toggle switch dynamics)"
    )

    TRIFLUOROMETHYL = FunctionalGroup(
        name="CF3",
        smarts="[C](F)(F)F",
        description="Trifluoromethyl (increases metabolic stability)"
    )

    PMB_ETHER = FunctionalGroup(
        name="PMB",
        smarts="[CH2]c1ccc(OC)cc1",
        description="para-Methoxybenzyl ether (CYP2D6-labile prodrug)"
    )

    PIVALOYL_ESTER = FunctionalGroup(
        name="Piv",
        smarts="C(=O)C(C)(C)C",
        description="Pivaloyl ester (carboxylesterase substrate)"
    )

    TRIAZOLE = FunctionalGroup(
        name="Triazole",
        smarts="c1nnn[nH]1",
        description="1,2,3-Triazole linker (click chemistry product)"
    )


@dataclass(frozen=True)
class PhysicochemicalProperties:
    """Immutable physicochemical property container."""
    molecular_weight: float  # g/mol
    logP: float  # Partition coefficient
    hbd: int  # H-bond donors
    hba: int  # H-bond acceptors
    tpsa: float  # Topological polar surface area (Å²)
    rotatable_bonds: int
    aromatic_rings: int

    def __post_init__(self):
        if self.molecular_weight <= 0:
            raise ValueError("Molecular weight must be positive")
        if self.hbd < 0 or self.hba < 0:
            raise ValueError("H-bond counts cannot be negative")

    @property
    def lipinski_violations(self) -> int:
        """Count Lipinski Rule of Five violations."""
        violations = 0
        if self.molecular_weight > 500:
            violations += 1
        if self.logP > 5:
            violations += 1
        if self.hbd > 5:
            violations += 1
        if self.hba > 10:
            violations += 1
        return violations

    @property
    def is_drug_like(self) -> bool:
        """Check if compound passes Lipinski's Rule of Five."""
        return self.lipinski_violations <= 1

    @property
    def bbb_permeable(self) -> bool:
        """Estimate blood-brain barrier permeability (Egan rule)."""
        return self.tpsa <= 90 and self.logP >= 1 and self.logP <= 5


@dataclass(frozen=True)
class BindingAffinity:
    """Type-safe binding affinity representation."""
    target: ReceptorTarget
    ki_nm: float  # Inhibition constant in nanomolar
    selectivity_ratio: Optional[float] = None  # vs primary off-target

    def __post_init__(self):
        if self.ki_nm <= 0:
            raise ValueError("Ki must be positive")

    @property
    def affinity_class(self) -> Literal["high", "moderate", "low"]:
        """Classify binding affinity."""
        if self.ki_nm < 10:
            return "high"
        elif self.ki_nm < 100:
            return "moderate"
        return "low"


@dataclass(frozen=True)
class BiasProfile:
    """Signaling bias quantification."""
    pathway: SignalingPathway
    efficacy: float  # 0-1 scale
    bias_factor: float  # Log bias vs reference

    @property
    def is_biased(self) -> bool:
        """Check if compound shows significant bias (>10-fold)."""
        return abs(self.bias_factor) >= 1.0  # log(10) = 1


@dataclass
class MolecularStructure(ABC):
    """Abstract base class for molecular structures."""
    smiles: SMILES
    name: str

    @abstractmethod
    def get_core_scaffold(self) -> SMILES:
        """Return the core scaffold SMILES."""
        pass


@dataclass
class THCScaffold(MolecularStructure):
    """Delta-9-THC scaffold representation."""
    c1_substituent: FunctionalGroup
    c3_chain: str
    c9_substituent: FunctionalGroup

    def get_core_scaffold(self) -> SMILES:
        """Return dibenzopyran core."""
        return SMILES("C1CC2C(C)(C)Oc3c(ccc(c3)-c3ccccc3)C2=CC1")


@dataclass
class ProdrugDesign:
    """Type-safe prodrug specification."""
    parent_smiles: SMILES
    prodrug_smiles: SMILES
    protecting_group: FunctionalGroup
    activation_enzyme: ActivationMechanism
    target_tissue: Literal["cortex", "hypothalamus", "systemic"]
    activation_rate: Optional[float] = None  # kcat in s^-1

    @property
    def is_cns_targeted(self) -> bool:
        """Check if prodrug targets CNS activation."""
        return self.target_tissue in ("cortex", "hypothalamus")


@dataclass
class MultiTargetProfile:
    """Multi-target ligand pharmacology."""
    primary_target: BindingAffinity
    secondary_targets: list[BindingAffinity]
    therapeutic_ratio: Optional[float] = None

    @property
    def is_selective(self) -> bool:
        """Check for >10-fold selectivity over off-targets."""
        if not self.secondary_targets:
            return True
        min_secondary_ki = min(t.ki_nm for t in self.secondary_targets)
        return self.primary_target.ki_nm / min_secondary_ki >= 10


@dataclass
class CompoundSpecification:
    """Complete type-safe compound specification."""
    compound_id: str
    strategy: CompoundStrategy
    structure: MolecularStructure
    physicochemical: PhysicochemicalProperties
    binding: list[BindingAffinity]
    signaling_bias: Optional[list[BiasProfile]] = None
    prodrug: Optional[ProdrugDesign] = None
    multi_target: Optional[MultiTargetProfile] = None
    modifications: list[str] = field(default_factory=list)
    rationale: list[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate compound specification."""
        if self.strategy == CompoundStrategy.PRODRUG_ACTIVATION and not self.prodrug:
            raise ValueError("Prodrug strategy requires ProdrugDesign")
        if self.strategy == CompoundStrategy.MULTI_TARGET and not self.multi_target:
            raise ValueError("Multi-target strategy requires MultiTargetProfile")

    @property
    def is_viable_lead(self) -> bool:
        """Check if compound meets lead viability criteria."""
        if not self.physicochemical.is_drug_like:
            return False
        if not any(b.ki_nm < 100 for b in self.binding if b.target == ReceptorTarget.CB1):
            return False
        return True


# Type aliases for common patterns
ModificationList = list[str]
RationaleList = list[str]
BindingProfile = list[BindingAffinity]
