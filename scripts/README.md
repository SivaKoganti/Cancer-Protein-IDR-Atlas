# Molecular Docking Workflow and Figure Generation Scripts

Reproducible computational and biophysical analysis pipeline for lead-HSA-5-fluorouracil binding study.

## Overview

This directory contains Python scripts for:
1. **Molecular docking simulations** (AutoDock Vina, GOLD)
2. **Allosteric pathway analysis**
3. **Publication-quality figure generation**
4. **Results summary reporting**

## Requirements

### System Requirements
- Python 3.8+
- ~2-4 GB disk space for docking results
- Molecular docking software (optional):
  - AutoDock Vina (free) - Required for docking
  - GOLD (CCDC) - Optional, requires license
  - OpenBabel - For ligand format conversion

### Python Packages

```bash
pip install numpy pandas matplotlib seaborn biopython scipy
```

### Optional (for full functionality)

```bash
pip install rdkit      # For ligand structure generation from SMILES
pip install meeko      # For PDBQT preparation (AutoDock-compatible)
```

## Installation

### 1. Install Molecular Docking Tools

**AutoDock Vina (Recommended - Free):**
```bash
# Ubuntu/Debian
sudo apt-get install autodock-vina

# macOS (Homebrew)
brew install autodock-vina

# Or download from: http://vina.scripps.edu/
```

**GOLD (Optional - Commercial):**
Requires CCDC license. Instructions: https://www.ccdc.cam.ac.uk/support/

### 2. Clone or Download Scripts

```bash
cd /path/to/Cancer-Protein-IDR-Atlas/scripts
```

### 3. Install Python Dependencies

```bash
pip install -r requirements.txt
```

## Usage

### Running the Molecular Docking Workflow

```bash
python molecular_docking_workflow.py \
  --work-dir ./docking_results \
  --receptor 1E7I.pdb
```

**Arguments:**
- `--work-dir`: Output directory (default: `./docking_results`)
- `--receptor`: HSA PDB file (default: `1E7I.pdb` - downloaded automatically)

**Output:**
```
docking_results/
├── 1E7I_receptor.pdbqt         # Prepared HSA receptor
├── lead.pdbqt                  # Lead ion (Pb2+)
├── 5fu.pdbqt                   # 5-FU ligand
├── vina_lead_output.pdbqt      # Vina lead docking results
├── vina_5fu_output.pdbqt       # Vina 5-FU docking results
├── docking_results.json        # Parsed results
└── docking_summary_report.txt  # Summary report
```

### Generating Publication Figures

```bash
python generate_figures.py --output-dir ./figures
```

**Output:**
```
figures/
├── Figure_3.1_Lead_Binding_Poses.png
├── Figure_3.2_5FU_Docking_Comparison.png
├── Figure_3.3_Allosteric_Pathway.png
├── Figure_3.4_Viscometry.png
├── Figure_3.5_Fluorescence_Spectroscopy.png
├── Figure_3.6_CD_Spectroscopy.png
├── Figure_3.7_DSC_Analysis.png
└── Figure_3.8_FTIR_Analysis.png
```

All figures generated at 300 dpi, publication quality.

## Detailed Script Documentation

### `molecular_docking_workflow.py`

Main docking pipeline orchestrating:

1. **Receptor Preparation** (Step 1)
   - Downloads HSA from PDB (1E7I.pdb)
   - Converts PDB → PDBQT with Gasteiger charges
   - Output: `*_receptor.pdbqt`

2. **Ligand Preparation** (Step 2)
   - Creates Pb²⁺ ion file
   - Generates 5-FU from SMILES (requires RDKit)
   - Output: `*.pdbqt` files

3. **Docking Box Definition** (Step 3)
   - 126 Å³ grid to encompass entire HSA
   - Center at protein centroid
   - Spacing: 0.375 Å

4. **AutoDock Vina Docking** (Step 4)
   - Lead: 20 conformations, biased toward Cys-34
   - 5-FU: 20 conformations, both native and lead-bound protein
   - Output: Binding free energies + poses

5. **Allosteric Pathway Analysis** (Step 5)
   - Identifies bridging residues: Cys-34 → Lys-129 → Asp-183 → Trp-214 → Lys-199
   - Calculates inter-residue distances
   - Computes contact frequencies

6. **Summary Report** (Step 6)
   - Binding energies table
   - Pathway analysis summary
   - Output: `docking_summary_report.txt`

### `generate_figures.py`

Publication-quality figure generator:

**Figure 3.1: Lead Binding Poses**
- Panel A: 4-coordinate Pb-Cys34 geometry
- Panel B: Ensemble of 20 docked poses (RMSD clustering)
- Panel C: Primary vs. secondary lead-binding sites

**Figure 3.2: 5-FU Docking Comparison**
- Panel A: Native HSA binding site distribution
- Panel B: Lead-bound HSA binding site distribution
- Panel C: Predicted binding free energy distributions
- Panel D: Conformational heterogeneity (RMSD analysis)

**Figure 3.3: Allosteric Pathway**
- Panel A: Inter-residue distances along pathway
- Panel B: Contact frequency heatmap (residue pairs)
- Panel C: Network visualization of allosteric pathway

**Figure 3.4: Viscometry**
- Specific viscosity (ηSP) across all treatment conditions
- Error bars from n=3 replicates
- Interpretation of conformational changes

**Figure 3.5: Fluorescence Spectroscopy**
- Panel A: Tryptophan fluorescence peak area quenching
- Panel B: Peak position shift (environmental hydrophilicity)
- Panel C: Stern-Volmer binding constant extraction
- Panel D: Stern-Volmer plots (native vs. lead-bound)

**Figure 3.6: CD Spectroscopy**
- Panel A: Far-UV CD spectra (190-260 nm)
- Panel B: Dose-dependent helical content loss

**Figure 3.7: DSC Analysis**
- Panel A: Thermal unfolding thermograms (25-40°C)
- Panel B: Enthalpy changes (ΔH) across conditions

**Figure 3.8: FTIR Spectroscopy**
- Panel A: Full IR spectra (4000-450 cm⁻¹)
- Panel B: Amide I region (1600-1690 cm⁻¹)
- Panel C: Amide III region (1229-1301 cm⁻¹)
- Panel D: Percent change from control

## Key Metrics & Expected Results

### Docking Results
| Condition | Binding Energy (kcal/mol) | Occupancy |
|---|---|---|
| Lead-HSA (Cys-34) | -7.8 ± 0.3 | 78% |
| 5-FU-HSA (native) | -6.3 ± 0.4 | 85% (site II) |
| 5-FU-HSA (lead-bound) | -4.5 ± 0.6 | 42% (site II) |
| **Affinity change (ΔΔG)** | **-1.8 kcal/mol** | **3-4 fold reduction** |

### Biophysical Validation
| Method | Finding | Interpretation |
|---|---|---|
| Viscometry | ηSP increases 2.1-3.4-fold | Protein elongation |
| Fluorescence | Binding constant reduction 2-4-fold | Reduced 5-FU affinity |
| CD Spectroscopy | Helix loss 10-31% | Secondary structure disruption |
| DSC | ΔH increases 2.8-fold | Compensatory tertiary packing |
| DR-FTIR | Amide I/III %R decrease 3.6-4.6% | Altered hydrogen bonding |

## Extending the Workflow

### Adding Other Metal Ions

Modify `prepare_ligands()` in `molecular_docking_workflow.py`:

```python
ligands = {
    "lead": "Pb2+",
    "cadmium": "Cd2+",
    "mercury": "Hg2+",
    "5fu": "CC(=O)NC1=C(C(=O)NC(=C1F)N)F"
}
```

### Adding Other HSA-Binding Drugs

```python
ligands = {
    "5fu": "CC(=O)NC1=C(C(=O)NC(=C1F)N)F",
    "warfarin": "CC(=O)CC(C1=CC=CC=C1)C2=C(O)C3=C(C=C2)C(=CC=C3)C(C)C",
    "ibuprofen": "CC(C)CC1=CC(=C(C=C1)C(C)C(=O)O)"
}
```

### Using Alternative Docking Engines

**GOLD Integration:**
```python
workflow.run_gold_docking(ligand_pdb, ligand_name, gold_dir="/opt/ccdc/gold")
```

**Requires:** CCDC GOLD license

## Troubleshooting

### Error: "meeko_prepare_receptor.py: command not found"
**Solution:** Install meeko package
```bash
pip install meeko
```

### Error: "vina: command not found"
**Solution:** Install AutoDock Vina
```bash
# Ubuntu/Debian
sudo apt-get install autodock-vina

# macOS
brew install autodock-vina
```

### Error: "ModuleNotFoundError: No module named 'rdkit'"
**Solution:** Install RDKit (optional, for 5-FU structure generation)
```bash
conda install -c conda-forge rdkit
# OR
pip install rdkit-pypi
```

### Docking produces no results
**Check:**
1. Receptor PDBQT file exists and is valid
2. Ligand PDBQT file exists and is valid
3. Vina installation: `which vina`
4. File permissions in working directory

## Performance & Runtime

| Step | Time | Hardware |
|---|---|---|
| Receptor preparation | ~2 min | Single core |
| Lead docking (20 modes) | ~5 min | Single core |
| 5-FU docking (20 modes) | ~10 min | Single core |
| Pathway analysis | <1 min | Single core |
| Figure generation | ~2 min | Single core |
| **Total** | **~20 min** | **Single core** |

## Output Interpretation

### `docking_summary_report.txt`

Example output:
```
================================================================================
MOLECULAR DOCKING WORKFLOW SUMMARY REPORT
Lead-HSA and 5-Fluorouracil Binding Study
================================================================================

Working Directory: ./docking_results
Receptor (HSA): 1E7I.pdb

LEAD DOCKING RESULTS (AutoDock Vina)
----------------------------------------
Number of conformations: 20
Best binding energy: -7.80 kcal/mol
Mean binding energy: -7.62 kcal/mol
Binding energies: ['-7.80', '-7.75', '-7.70', '-7.68', '-7.65']...

5-FLUOROURACIL DOCKING RESULTS (AutoDock Vina)
----------------------------------------
Number of conformations: 20
Best binding energy: -6.30 kcal/mol
Mean binding energy: -6.18 kcal/mol
Binding energies: ['-6.30', '-6.25', '-6.20', '-6.15', '-6.10']...

ALLOSTERIC PATHWAY ANALYSIS
----------------------------------------
Pathway: Cys-34 (Pb) → Lys-199 (5-FU)
Total distance: 30.4 Å
Key residues: Cys-34, Lys-129, Asp-183, Trp-214, Lys-199
```

## Citation

If using these scripts, please cite:
```
Sivan et al. "Molecular Docking and Biophysical Studies of Lead Interference 
with Human Serum Albumin and 5-Fluorouracil Binding" [Journal], [Year].

Code: https://github.com/SivaKoganti/Cancer-Protein-IDR-Atlas
```

## Contributing

To extend these scripts:
1. Fork the repository
2. Create a feature branch
3. Add your modifications with documentation
4. Submit a pull request

## Support & Issues

For bugs, feature requests, or questions:
- GitHub Issues: https://github.com/SivaKoganti/Cancer-Protein-IDR-Atlas/issues
- Email: [contact information]

## License

MIT License - See LICENSE file for details

## Acknowledgments

- AutoDock Vina: Morris et al.
- GOLD: CCDC
- RDKit: For cheminformatics
- Matplotlib, Seaborn: For visualization

---

**Last Updated:** 2026-08-20
**Status:** Ready for production use
