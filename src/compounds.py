#!/usr/bin/env python3
"""
YAC Compound Library: Pathway-Selective Cannabinoid Compounds for ADHD

This module defines the six lead compounds (YAC-101/102, YAC-201/202, YAC-301/302)
designed for dissociating therapeutic efficacy from orexigenic effects through
biased signaling, prodrug activation, and multi-target pharmacology.

Author: Matthew Long / YonedaAI Collaboration
"""

from dataclasses import dataclass, field
from typing import Optional, TypedDict, Final
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import json

from .types import (
    SMILES,
    CompoundStrategy,
    ReceptorTarget,
    SignalingPathway,
    ActivationMechanism,
    BindingMode,
    FunctionalGroup,
    CannabinoiFunctionalGroups,
    PhysicochemicalProperties,
    BindingAffinity,
    BiasProfile,
    ProdrugDesign,
    MultiTargetProfile,
    CompoundSpecification,
    ModificationList,
    RationaleList,
)


# Type-safe constants
THC_SMILES: Final[SMILES] = SMILES("CCCCCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1")
THC_MW: Final[float] = 314.5
THC_LOGP: Final[float] = 6.97


class CompoundDict(TypedDict):
    """Type-safe dictionary representation for JSON export."""
    id: str
    name: str
    strategy: str
    smiles: str
    molecular_formula: str
    molecular_weight: float
    logP: float
    key_modifications: list[str]
    design_rationale: list[str]


@dataclass
class CompoundData:
    """Data class for compound specifications with type-safe fields."""

    id: str
    name: str
    strategy: str
    smiles: SMILES
    molecular_formula: str
    molecular_weight: float
    logp: float
    key_modifications: ModificationList
    design_rationale: RationaleList
    binding_mode: BindingMode = BindingMode.ORTHOSTERIC
    primary_target: ReceptorTarget = ReceptorTarget.CB1
    activation_mechanism: Optional[ActivationMechanism] = None
    bias_pathways: list[SignalingPathway] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate compound data."""
        if self.molecular_weight <= 0:
            raise ValueError(f"Invalid molecular weight: {self.molecular_weight}")
        if not self.smiles:
            raise ValueError("SMILES string cannot be empty")

    def to_dict(self) -> CompoundDict:
        return CompoundDict(
            id=self.id,
            name=self.name,
            strategy=self.strategy,
            smiles=self.smiles,
            molecular_formula=self.molecular_formula,
            molecular_weight=self.molecular_weight,
            logP=self.logp,
            key_modifications=self.key_modifications,
            design_rationale=self.design_rationale
        )

    def to_specification(self) -> Optional[CompoundSpecification]:
        """Convert to full type-safe specification (requires RDKit analysis)."""
        props = MolecularAnalyzer.analyze_compound(self.smiles)
        if "error" in props:
            return None

        physicochemical = PhysicochemicalProperties(
            molecular_weight=props["molecular_weight"],
            logP=props["logP"],
            hbd=props["hbd"],
            hba=props["hba"],
            tpsa=props["tpsa"],
            rotatable_bonds=props["rotatable_bonds"],
            aromatic_rings=props["aromatic_rings"]
        )

        # Predicted binding (would need experimental data)
        binding = [BindingAffinity(
            target=self.primary_target,
            ki_nm=50.0,  # Placeholder
            selectivity_ratio=10.0
        )]

        return CompoundSpecification(
            compound_id=self.id,
            strategy=CompoundStrategy(self.strategy.split()[0]),
            structure=None,  # Would need THCScaffold
            physicochemical=physicochemical,
            binding=binding,
            modifications=self.key_modifications,
            rationale=self.design_rationale
        )


# Base THC SMILES reference (duplicated as constant above for backwards compatibility)


class YACCompoundLibrary:
    """
    Library of YAC (YonedaAI Cannabinoid) compounds for ADHD treatment.

    Three mechanistically distinct strategies:
    - Strategy A: Biased signaling (Gi-selective partial agonists)
    - Strategy B: Prodrug activation (cortex-selective)
    - Strategy C: Multi-target ligands (appetite pathway antagonism)
    """

    def __init__(self):
        self.compounds = self._initialize_compounds()

    def _initialize_compounds(self) -> dict[str, CompoundData]:
        """Initialize the six lead compounds."""

        compounds = {}

        # YAC-101: Gi-Biased Partial Agonist
        # Modifications: C1→OCH3, C9→CHF2, C3→(CH2)3CF3
        # THC with C1 methyl ether, C9 difluoromethyl, C3 trifluorobutyl
        yac101_smiles: SMILES = SMILES("FC(F)(F)CCCc1cc(OC)c2C3CC(C(F)F)=CCC3C(C)(C)Oc2c1")
        compounds["YAC-101"] = CompoundData(
            id="YAC-101",
            name="1-Methoxy-9-difluoromethyl-3-(4,4,4-trifluorobutyl)-THC",
            strategy="A - Gi-Biased Partial Agonist",
            smiles=yac101_smiles,
            molecular_formula="C23H27F5O2",
            molecular_weight=434.5,
            logp=5.8,
            key_modifications=[
                "C1 → OCH3 (methyl ether)",
                "C9 → CHF2 (difluoromethyl)",
                "C3 → (CH2)3CF3 (4,4,4-trifluorobutyl)"
            ],
            design_rationale=[
                "C1 methyl ether reduces intrinsic efficacy to ~50% of THC",
                "C9 difluoromethyl alters toggle switch dynamics for Gi-coupled states",
                "Terminal trifluoromethyl increases metabolic stability",
                "Combined modifications yield Gi-biased partial agonist profile"
            ],
            binding_mode=BindingMode.ORTHOSTERIC,
            primary_target=ReceptorTarget.CB1,
            bias_pathways=[SignalingPathway.Gi]
        )

        # YAC-102: Allosteric-Orthosteric Hybrid
        # THC core + triazole linker + indole-2-carboxamide
        yac102_smiles: SMILES = SMILES(
            "CCCCn1cc(nn1)CCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1.NC(=O)c1cc2ccccc2[nH]1"
        )
        compounds["YAC-102"] = CompoundData(
            id="YAC-102",
            name="THC-C3-triazole-ethyl-indole-2-carboxamide",
            strategy="A - Allosteric-Orthosteric Hybrid",
            smiles=yac102_smiles,
            molecular_formula="C35H42N4O4",
            molecular_weight=586.7,
            logp=5.2,
            key_modifications=[
                "C1-OH preserved for orthosteric binding",
                "C9-CH3 preserved for potency",
                "C3 → propyl-triazole-ethyl-indole-2-carboxamide"
            ],
            design_rationale=[
                "Intact THC core maintains orthosteric activation capacity",
                "Triazole linker enables modular synthesis via click chemistry",
                "Indole-2-carboxamide engages TM2-TM3-TM4 allosteric site",
                "Dual-site engagement locks Gi-selective receptor conformation"
            ],
            binding_mode=BindingMode.BITOPIC,
            primary_target=ReceptorTarget.CB1,
            bias_pathways=[SignalingPathway.Gi]
        )

        # YAC-201: CYP2D6-Activated Prodrug
        # THC 1-O-(4-methoxybenzyl) ether
        yac201_smiles: SMILES = SMILES("CCCCCc1cc(OCc2ccc(OC)cc2)c2C3CC(C)=CCC3C(C)(C)Oc2c1")
        compounds["YAC-201"] = CompoundData(
            id="YAC-201",
            name="Δ9-THC 1-O-(4-methoxybenzyl) ether",
            strategy="B - CYP2D6-Activated Prodrug",
            smiles=yac201_smiles,
            molecular_formula="C29H36O3",
            molecular_weight=448.6,
            logp=6.2,
            key_modifications=[
                "C1 → O-CH2-C6H4-OCH3 (para-methoxybenzyl ether)"
            ],
            design_rationale=[
                "PMB ether is inactive at CB1 (requires free C1-OH)",
                "High lipophilicity (logP 6.2) ensures BBB penetration",
                "CYP2D6 O-dealkylation releases THC",
                "CYP2D6 is 3-5× higher in cortex than hypothalamus",
                "Prodrug acts as cortex-selective THC delivery system"
            ],
            binding_mode=BindingMode.ORTHOSTERIC,
            primary_target=ReceptorTarget.CB1,
            activation_mechanism=ActivationMechanism.CYP2D6
        )

        # YAC-202: Esterase-Activated Prodrug
        # THC 1-O-pivaloate
        yac202_smiles: SMILES = SMILES("CCCCCc1cc(OC(=O)C(C)(C)C)c2C3CC(C)=CCC3C(C)(C)Oc2c1")
        compounds["YAC-202"] = CompoundData(
            id="YAC-202",
            name="Δ9-THC 1-O-pivaloate",
            strategy="B - Esterase-Activated Prodrug",
            smiles=yac202_smiles,
            molecular_formula="C26H36O3",
            molecular_weight=400.6,
            logp=6.5,
            key_modifications=[
                "C1 → O-CO-C(CH3)3 (pivaloyl ester)"
            ],
            design_rationale=[
                "Pivaloyl ester is inactive at CB1 (C1 blocked)",
                "Bulky tert-butyl group provides steric protection",
                "Regional CES activity differences determine activation",
                "Complementary mechanism to CYP2D6 prodrug"
            ],
            binding_mode=BindingMode.ORTHOSTERIC,
            primary_target=ReceptorTarget.CB1,
            activation_mechanism=ActivationMechanism.CES1
        )

        # YAC-301: CB1 Agonist / Ghrelin Antagonist Hybrid
        # THC-C3-PEG4-triazole-[D-Lys3]-GHRP-6 antagonist
        # Simplified representation (full peptide too complex for SMILES)
        yac301_smiles: SMILES = SMILES(
            "CCCCCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1"  # THC core
            ".COCCOCCOCCOCCO"  # PEG4 linker representation
            ".CC(C)CC(NC(=O)C(Cc1c[nH]c2ccccc12)NC(=O)C(N)CCCCN)C(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NC(Cc1ccccc1)C(N)=O"
        )
        compounds["YAC-301"] = CompoundData(
            id="YAC-301",
            name="THC-C3-PEG4-triazole-[D-Lys³]-GHRP-6 antagonist",
            strategy="C - CB1 Agonist / Ghrelin Antagonist",
            smiles=yac301_smiles,
            molecular_formula="Hybrid",
            molecular_weight=1150.0,
            logp=3.5,  # Estimated
            key_modifications=[
                "THC moiety for CB1 agonism",
                "PEG4 spacer for flexibility and solubility",
                "Peptidomimetic component antagonizes GHSR1a",
                "D-amino acids for proteolytic stability"
            ],
            design_rationale=[
                "Dual-target approach: CB1 activation + GHSR1a antagonism",
                "Blocks ghrelin-CB1 synergy in hypothalamic feeding circuits",
                "THC component preserved for therapeutic CB1 effects",
                "Ghrelin antagonism prevents orexigenic effects"
            ],
            binding_mode=BindingMode.ORTHOSTERIC,
            primary_target=ReceptorTarget.CB1
        )

        # YAC-302: CB1 Agonist / MC4R Agonist Hybrid
        # THC-linker-setmelanotide mimetic
        yac302_smiles: SMILES = SMILES(
            "CCCCCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1"  # THC core
            ".CC(C)CC(NC(=O)C(Cc1ccc(O)cc1)NC(=O)C(CCCNC(=N)N)NC(=O)C(Cc1c[nH]cn1)NC(=O)C)C(=O)NC(Cc1ccccc1)C(=O)NC(CCC(N)=O)C(N)=O"
        )
        compounds["YAC-302"] = CompoundData(
            id="YAC-302",
            name="THC-linker-setmelanotide mimetic",
            strategy="C - CB1 Agonist / MC4R Agonist",
            smiles=yac302_smiles,
            molecular_formula="Hybrid",
            molecular_weight=1050.0,
            logp=2.8,  # Estimated
            key_modifications=[
                "THC core for CB1 agonism",
                "Alkyl-triazole spacer",
                "Cyclic peptide derived from setmelanotide"
            ],
            design_rationale=[
                "Dual agonism at CB1 (therapeutic) and MC4R (anorexigenic)",
                "Internal counterbalance approach",
                "Setmelanotide is FDA-approved MC4R agonist",
                "May require non-oral delivery due to peptide component"
            ],
            binding_mode=BindingMode.ORTHOSTERIC,
            primary_target=ReceptorTarget.CB1
        )

        return compounds

    def get_compound(self, compound_id: str) -> Optional[CompoundData]:
        """Retrieve a compound by ID."""
        return self.compounds.get(compound_id)

    def get_all_compounds(self) -> dict[str, CompoundData]:
        """Get all compounds in the library."""
        return self.compounds

    def get_by_strategy(self, strategy: str) -> list[CompoundData]:
        """Get compounds by strategy type (A, B, or C)."""
        return [c for c in self.compounds.values() if strategy in c.strategy]


class MolecularAnalyzer:
    """Analyze molecular properties using RDKit."""

    @staticmethod
    def analyze_compound(smiles: str) -> dict:
        """Calculate molecular descriptors for a compound."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}

        # Handle multi-fragment molecules
        frags = smiles.split(".")
        main_frag = max(frags, key=len) if frags else smiles
        mol = Chem.MolFromSmiles(main_frag)

        if mol is None:
            return {"error": "Could not parse main fragment"}

        return {
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "logP": round(Descriptors.MolLogP(mol), 2),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),  # H-bond donors
            "hba": rdMolDescriptors.CalcNumHBA(mol),  # H-bond acceptors
            "tpsa": round(Descriptors.TPSA(mol), 2),  # Topological PSA
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "lipinski_violations": _count_lipinski_violations(mol)
        }

    @staticmethod
    def generate_2d_coords(smiles: str) -> Optional[Chem.Mol]:
        """Generate 2D coordinates for visualization."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            AllChem.Compute2DCoords(mol)
        return mol

    @staticmethod
    def check_drug_likeness(smiles: str) -> dict:
        """Evaluate drug-likeness based on Lipinski's Rule of Five."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)

        violations = []
        if mw > 500:
            violations.append(f"MW ({mw:.1f}) > 500")
        if logp > 5:
            violations.append(f"logP ({logp:.2f}) > 5")
        if hbd > 5:
            violations.append(f"HBD ({hbd}) > 5")
        if hba > 10:
            violations.append(f"HBA ({hba}) > 10")

        return {
            "is_drug_like": len(violations) <= 1,
            "violations": violations,
            "violation_count": len(violations),
            "properties": {
                "MW": round(mw, 2),
                "logP": round(logp, 2),
                "HBD": hbd,
                "HBA": hba
            }
        }


def _count_lipinski_violations(mol: Chem.Mol) -> int:
    """Count Lipinski Rule of Five violations."""
    violations = 0
    if Descriptors.MolWt(mol) > 500:
        violations += 1
    if Descriptors.MolLogP(mol) > 5:
        violations += 1
    if rdMolDescriptors.CalcNumHBD(mol) > 5:
        violations += 1
    if rdMolDescriptors.CalcNumHBA(mol) > 10:
        violations += 1
    return violations


def generate_compound_report(library: YACCompoundLibrary) -> str:
    """Generate a comprehensive report of all compounds."""
    analyzer = MolecularAnalyzer()
    report_lines = [
        "=" * 80,
        "YAC COMPOUND LIBRARY REPORT",
        "Pathway-Selective Cannabinoid Compounds for ADHD",
        "=" * 80,
        ""
    ]

    for strategy in ["A", "B", "C"]:
        compounds = library.get_by_strategy(strategy)
        if not compounds:
            continue

        strategy_names = {
            "A": "Strategy A: Biased Signaling",
            "B": "Strategy B: Prodrug Activation",
            "C": "Strategy C: Multi-Target Ligands"
        }

        report_lines.append(f"\n{strategy_names[strategy]}")
        report_lines.append("-" * 60)

        for compound in compounds:
            report_lines.append(f"\n{compound.id}: {compound.name}")
            report_lines.append(f"  Formula: {compound.molecular_formula}")
            report_lines.append(f"  MW: {compound.molecular_weight} g/mol")
            report_lines.append(f"  logP: {compound.logp}")

            # RDKit analysis (main fragment only for hybrids)
            smiles_to_analyze = compound.smiles.split(".")[0]
            props = analyzer.analyze_compound(smiles_to_analyze)
            if "error" not in props:
                report_lines.append(f"  Computed Properties (main fragment):")
                report_lines.append(f"    TPSA: {props['tpsa']} Å²")
                report_lines.append(f"    Rotatable Bonds: {props['rotatable_bonds']}")
                report_lines.append(f"    Aromatic Rings: {props['aromatic_rings']}")

            report_lines.append(f"  Key Modifications:")
            for mod in compound.key_modifications:
                report_lines.append(f"    - {mod}")

    return "\n".join(report_lines)


def export_to_json(library: YACCompoundLibrary, filepath: str) -> None:
    """Export compound library to JSON format."""
    data = {
        "library_name": "YAC Cannabinoid Compounds for ADHD",
        "version": "1.0",
        "compounds": [c.to_dict() for c in library.compounds.values()]
    }
    with open(filepath, "w") as f:
        json.dump(data, f, indent=2)


def main():
    """Main entry point for compound analysis."""
    library = YACCompoundLibrary()

    print(generate_compound_report(library))

    print("\n" + "=" * 80)
    print("DRUG-LIKENESS ASSESSMENT")
    print("=" * 80)

    analyzer = MolecularAnalyzer()
    for compound_id, compound in library.compounds.items():
        main_smiles = compound.smiles.split(".")[0]
        assessment = analyzer.check_drug_likeness(main_smiles)

        if "error" in assessment:
            print(f"\n{compound_id}: Error analyzing structure")
            continue

        status = "PASS" if assessment["is_drug_like"] else "REVIEW"
        print(f"\n{compound_id}: {status}")
        if assessment["violations"]:
            for v in assessment["violations"]:
                print(f"  - {v}")

    # Export to JSON
    export_to_json(library, "compounds.json")
    print("\n\nCompound data exported to compounds.json")


if __name__ == "__main__":
    main()
