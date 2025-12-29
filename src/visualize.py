#!/usr/bin/env python3
"""
Visualization module for YAC Compound Library.

Generates molecular structure images and comparison visualizations.
"""

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from compounds import YACCompoundLibrary


def draw_compound_grid(library: YACCompoundLibrary, output_path: str) -> None:
    """Draw all compounds in a grid layout."""
    compounds = list(library.compounds.values())
    mols = []
    legends = []

    for compound in compounds:
        # Use main fragment for multi-fragment molecules
        smiles = compound.smiles.split(".")[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
            legends.append(compound.id)

    if mols:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=3,
            subImgSize=(400, 400),
            legends=legends
        )
        img.save(output_path)
        print(f"Compound grid saved to {output_path}")


def draw_single_compound(smiles: str, compound_id: str, output_path: str) -> None:
    """Draw a single compound structure."""
    main_smiles = smiles.split(".")[0]
    mol = Chem.MolFromSmiles(main_smiles)

    if mol is None:
        print(f"Error: Could not parse SMILES for {compound_id}")
        return

    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    with open(output_path, "wb") as f:
        f.write(drawer.GetDrawingText())
    print(f"Structure for {compound_id} saved to {output_path}")


def draw_thc_comparison(library: YACCompoundLibrary, output_path: str) -> None:
    """Draw THC alongside Strategy B prodrugs for comparison."""
    thc_smiles = "CCCCCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1"
    thc_mol = Chem.MolFromSmiles(thc_smiles)

    yac201 = library.get_compound("YAC-201")
    yac202 = library.get_compound("YAC-202")

    mols = [thc_mol]
    legends = ["Δ9-THC (parent)"]

    if yac201:
        mol = Chem.MolFromSmiles(yac201.smiles.split(".")[0])
        if mol:
            mols.append(mol)
            legends.append("YAC-201 (CYP2D6 prodrug)")

    if yac202:
        mol = Chem.MolFromSmiles(yac202.smiles.split(".")[0])
        if mol:
            mols.append(mol)
            legends.append("YAC-202 (esterase prodrug)")

    for mol in mols:
        AllChem.Compute2DCoords(mol)

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=3,
        subImgSize=(400, 400),
        legends=legends
    )
    img.save(output_path)
    print(f"THC comparison saved to {output_path}")


def generate_all_visualizations(output_dir: str = "images") -> None:
    """Generate all visualization outputs."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    library = YACCompoundLibrary()

    # Grid of all compounds
    draw_compound_grid(library, str(output_path / "yac_compound_grid.png"))

    # Individual compound structures
    for compound_id, compound in library.compounds.items():
        draw_single_compound(
            compound.smiles,
            compound_id,
            str(output_path / f"{compound_id.lower()}_structure.png")
        )

    # THC vs prodrugs comparison
    draw_thc_comparison(library, str(output_path / "thc_prodrug_comparison.png"))

    print(f"\nAll visualizations saved to {output_dir}/")


if __name__ == "__main__":
    generate_all_visualizations()
