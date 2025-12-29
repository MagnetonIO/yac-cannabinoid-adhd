# YAC Cannabinoid Compound Library for ADHD

Rational Design of Pathway-Selective Cannabinoid Compounds for Attention-Deficit/Hyperactivity Disorder: Dissociating Therapeutic Efficacy from Orexigenic Effects Through Biased Signaling, Prodrug Activation, and Multi-Target Pharmacology.

## Overview

This repository contains the computational chemistry implementation and research paper for six lead compounds (YAC-101/102, YAC-201/202, YAC-301/302) designed to preserve ADHD-relevant prefrontal cortex modulation while eliminating appetite stimulation through three mechanistically distinct strategies:

| Strategy | Approach | Compounds |
|----------|----------|-----------|
| **A** | Biased Signaling | YAC-101, YAC-102 |
| **B** | Prodrug Activation | YAC-201, YAC-202 |
| **C** | Multi-Target Ligands | YAC-301, YAC-302 |

## Compound Library

### Strategy A: Biased Signaling

- **YAC-101**: Gi-Biased Partial Agonist (MW: 434.5, logP: 5.8)
- **YAC-102**: Allosteric-Orthosteric Hybrid (MW: 586.7, logP: 5.2)

### Strategy B: Prodrug Activation

- **YAC-201**: CYP2D6-Activated Prodrug (MW: 448.6, logP: 6.2)
- **YAC-202**: Esterase-Activated Prodrug (MW: 400.6, logP: 6.5)

### Strategy C: Multi-Target Ligands

- **YAC-301**: CB1 Agonist / Ghrelin Antagonist (MW: ~1150)
- **YAC-302**: CB1 Agonist / MC4R Agonist (MW: ~1050)

## Installation

```bash
# Clone the repository
git clone https://github.com/mlong/yac-cannabinoid-adhd.git
cd yac-cannabinoid-adhd

# Install dependencies
pip install -r requirements.txt

# Or install as package
pip install -e .
```

## Usage

### Analyze Compounds

```python
from src.compounds import YACCompoundLibrary, MolecularAnalyzer

# Initialize library
library = YACCompoundLibrary()

# Get specific compound
yac101 = library.get_compound("YAC-101")
print(f"Name: {yac101.name}")
print(f"SMILES: {yac101.smiles}")

# Analyze molecular properties
analyzer = MolecularAnalyzer()
props = analyzer.analyze_compound(yac101.smiles)
print(f"Computed MW: {props['molecular_weight']}")
print(f"TPSA: {props['tpsa']} Å²")

# Drug-likeness assessment
assessment = analyzer.check_drug_likeness(yac101.smiles)
print(f"Drug-like: {assessment['is_drug_like']}")
```

### Generate Visualizations

```python
from src.visualize import generate_all_visualizations

# Generate all compound structure images
generate_all_visualizations(output_dir="images")
```

### Command Line

```bash
# Run compound analysis
python src/compounds.py
```

## Repository Structure

```
yac-cannabinoid-adhd/
├── paper/
│   └── pathway-selective-cannabinoids-adhd.tex
├── src/
│   ├── __init__.py
│   ├── compounds.py      # Core compound definitions
│   └── visualize.py      # Molecular visualization
├── docs/
│   └── index.html        # GitHub Pages
├── requirements.txt
├── pyproject.toml
└── README.md
```

## Research Paper

The full research paper is available in the `paper/` directory:
- LaTeX source: `paper/pathway-selective-cannabinoids-adhd.tex`
- PDF: `paper/pathway-selective-cannabinoids-adhd.pdf`

### Citation

```bibtex
@article{long2024cannabinoid,
  title={Rational Design of Pathway-Selective Cannabinoid Compounds for
         Attention-Deficit/Hyperactivity Disorder},
  author={Long, Matthew and {YonedaAI Collaboration}},
  year={2024},
  publisher={YonedaAI Research}
}
```

## Key Scientific Concepts

### Biased Signaling (Strategy A)
Exploits differential CB1 receptor pathway engagement. Gi-biased compounds preferentially activate therapeutic pathways while minimizing β-arrestin-mediated appetite signaling.

### Prodrug Activation (Strategy B)
Utilizes regional differences in CYP2D6 and carboxylesterase expression for cortex-selective drug activation. CYP2D6 is 3-5× higher in cortex than hypothalamus.

### Multi-Target Pharmacology (Strategy C)
Simultaneously engages CB1 for therapeutic effects while antagonizing appetite-related receptors (GHSR1a, NPY Y1) or agonizing anorexigenic receptors (MC4R).

## Dependencies

- Python 3.10+
- RDKit 2023.9.1+
- NumPy
- Pillow

## License

MIT License

## Authors

- Matthew Long ([@mlong](https://github.com/mlong)) - Independent Researcher
- The YonedaAI Collaboration

## Acknowledgments

This research was conducted as part of the YonedaAI Research Collective's drug discovery initiative.
