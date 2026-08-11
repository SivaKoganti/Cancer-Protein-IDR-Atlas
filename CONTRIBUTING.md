# ⚠️ CONFIDENTIAL - NOT ACCEPTING PUBLIC CONTRIBUTIONS

This project is currently under confidential development. **Public contributions are NOT accepted** during this phase.

---

## Authorized Team Members Only

If you are an **authorized co-author or collaborator**, contact the corresponding author:

📧 **Siva Koganti** - s.koganti@[institution].edu

---

## Public Contributions (After Publication)

Once this work is published and licensing is established, contribution guidelines will be made available in a public release.

**Expected Timeline:** Late 2026 or early 2027 (following journal publication)

---

## Previous Version

*The following contribution guidelines apply ONLY after public release:*

---

## Code of Conduct

We are committed to providing a welcoming and inspiring community for all. Please read and follow our Code of Conduct:

- **Be respectful:** Treat all contributors with courtesy
- **Be inclusive:** Welcome diverse perspectives and backgrounds
- **Be constructive:** Provide helpful feedback and suggestions
- **Report issues:** Contact maintainers confidentially for concerns

---

## Ways to Contribute

### 1. Report Bugs 🐛

Found an error in the data or code? Help us fix it!

**How to report:**
1. Open an issue on GitHub: https://github.com/SivaKoganti/Cancer-Protein-IDR-Atlas/issues
2. Include:
   - **Title:** Brief description (e.g., "TP53 position 100 IDR score incorrect")
   - **Description:** What went wrong and expected behavior
   - **Steps to reproduce:** Code/command to trigger the bug
   - **Environment:** Python version, OS, conda environment
   - **Screenshots:** If applicable (e.g., incorrect visualization)
3. Label as `bug`

**Example Issue:**
```
Title: KRAS_atlas_with_clinvar.tsv missing conservation scores

Description:
When I load the KRAS atlas file, the conservation column contains NaN values 
for positions 1-50, while other genes have valid scores.

Steps to reproduce:
```python
import pandas as pd
kras = pd.read_csv('results/atlas/KRAS_atlas_with_clinvar.tsv', sep='\t')
print(kras[kras['conservation'].isna()].shape)  # Returns (50, 11)
```

Environment:
- Python 3.12.1
- Ubuntu 24.04 LTS
- conda env: cancer-idr-atlas

Label: bug
```

### 2. Suggest Enhancements 🚀

Have an idea to improve the atlas?

**How to suggest:**
1. Open an issue with label `enhancement`
2. Describe:
   - **What:** New feature or improvement
   - **Why:** Problem it solves or value it adds
   - **How:** Proposed implementation (optional)
   - **Examples:** Use cases or references

**Example Enhancement:**
```
Title: Add AlphaFold2 structural confidence scores to atlas

Description:
Recent AlphaFold2 releases provide pLDDT confidence scores for all human 
proteins. These would complement our existing structural proxy scores and 
provide experimental validation targets.

Value:
- Distinguish high-confidence IDR predictions from low-confidence
- Prioritize regions for experimental validation
- Improve disease mechanism interpretation

Implementation:
- Download pLDDT scores from EBI AlphaFold database
- Merge into per-residue atlas at build_atlas stage
- Add plddt_score column to atlas TSVs
- Update dashboard to color-code by pLDDT confidence

References:
- https://alphafold.ebi.ac.uk/
- Varadi et al. 2021, Nucleic Acids Res.

Label: enhancement
```

### 3. Improve Documentation 📚

Help other researchers use the atlas by improving docs!

**What to improve:**
- Clarify confusing sections in README or Supplementary Methods
- Add troubleshooting guides
- Write tutorials for specific analyses
- Add code examples
- Fix typos

**How to submit:**
1. Fork the repository
2. Edit `.md` files directly
3. Test locally (ensure formatting renders correctly)
4. Submit pull request with description of improvements

### 4. Add New Genes or Data 🧬

Want to expand the atlas beyond 50 CGC genes?

**Process:**
1. **Propose:** Open an issue describing which genes to add and why
2. **Validate:** Ensure data quality (see Data Standards below)
3. **Implement:** 
   - Add gene symbols to `config.yaml`
   - Run: `snakemake --cores 4 --use-conda`
   - Verify outputs: `ls results/atlas/{GENE}_atlas_with_clinvar.tsv`
4. **Test:** Ensure TSV files are complete and no NaN values
5. **Submit:** Pull request with updated config and validation

**Data Standards:**
- UniProt sequence must exist (canonical, not isoform)
- Genomic locus must be in GRCh38
- ClinVar variants must be available
- All per-residue scores must be computed (no NaN)

### 5. Fix Code Issues and Performance 🔧

See a potential bug or performance improvement in the code?

**Process:**
1. Create an issue describing the problem and proposed fix
2. Fork the repository
3. Create a branch: `git checkout -b fix/issue-description`
4. Make changes with clear commit messages
5. Add/update tests
6. Submit pull request

**Code Style:**
- Follow PEP 8 (Python)
- Add docstrings to functions
- Include type hints where possible
- Comment complex logic
- Keep lines < 100 characters

**Example:**
```python
def compute_disorder_enrichment(atlas_df, variant_col='clinvar_significance'):
    """
    Compute odds ratio of pathogenic variants in disordered regions.
    
    Parameters
    ----------
    atlas_df : pd.DataFrame
        Per-residue atlas with columns: iupred_score, clinvar_significance
    variant_col : str
        Column name for variant significance classification
    
    Returns
    -------
    float
        Odds ratio (pathogenic in IDR / pathogenic in ordered)
    """
    # Implementation...
    return odds_ratio
```

### 6. Contribute Analysis Code 📊

Developed a novel analysis using the atlas data?

**How to share:**
1. Create a new file in `scripts/analysis/` with clear name: `scripts/analysis/analyze_variant_enrichment.py`
2. Include:
   - Docstring explaining analysis
   - Input data requirements
   - Output format
   - Example usage
3. Add entry to `scripts/README.md`
4. Submit pull request

**Example Template:**
```python
#!/usr/bin/env python3
"""
Analyze variant enrichment in disordered regions.

Usage:
    python3 scripts/analysis/variant_enrichment_analysis.py \
        --atlas results/atlas/ \
        --output results/analysis/variant_enrichment.csv

Output:
    CSV file with columns: gene, disorder_fraction, variant_count, 
                          pathogenic_in_disorder, odds_ratio
"""

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--atlas', required=True, help='Atlas directory')
    parser.add_argument('--output', required=True, help='Output file')
    args = parser.parse_args()
    
    # Implementation...
```

---

## Pull Request Process

### Before You Start

1. **Check for existing PRs:** Don't duplicate work—search open PRs first
2. **Create an issue first:** Discuss major changes before implementing
3. **Keep it focused:** Each PR should address one issue/feature
4. **Update tests:** Ensure existing tests still pass
5. **Update docs:** If changing behavior, update README/methods

### Submission Steps

1. **Fork** the repository
2. **Branch:** Create a descriptive branch name
   - Bug fix: `fix/short-description`
   - Feature: `feature/short-description`
   - Documentation: `docs/short-description`
3. **Commit:** Write clear, concise commit messages
   ```
   Fix: Correct CDKN2A isoform matching in build_disorder_atlas.py
   
   - Changed gene name normalization to handle _1, _2 isoforms
   - Updated group_by_gene() to use regex matching
   - Verified all 50 genes now produce atlas output
   - Added unit tests for isoform edge cases
   ```
4. **Push:** Push to your fork
5. **Pull Request:**
   - Fill out PR template completely
   - Reference related issues: "Fixes #42"
   - Describe changes and testing
   - Include before/after if applicable
6. **Review:** Respond to reviewer feedback promptly

### PR Checklist

- [ ] Branch created from `mvp-initial`
- [ ] Changes address the issue/feature described
- [ ] Code follows PEP 8 style guide
- [ ] Added/updated docstrings and comments
- [ ] Tests added or updated (run: `pytest`)
- [ ] No new warnings or deprecated code
- [ ] README/docs updated if needed
- [ ] Commit messages are clear and concise

### Merge Criteria

Your PR will be merged when:
1. ✓ All tests pass
2. ✓ Code review approved (at least 1 maintainer)
3. ✓ No merge conflicts
4. ✓ Documentation updated
5. ✓ Changes follow project standards

---

## Development Setup

### Local Development Environment

```bash
# 1. Clone your fork
git clone https://github.com/YOUR_USERNAME/Cancer-Protein-IDR-Atlas.git
cd Cancer-Protein-IDR-Atlas

# 2. Add upstream (original repo)
git remote add upstream https://github.com/SivaKoganti/Cancer-Protein-IDR-Atlas.git

# 3. Create conda environment
conda env create -f environment.yml
conda activate cancer-idr-atlas

# 4. Install in development mode
pip install -e .

# 5. Run tests
pytest -v

# 6. Create feature branch
git checkout -b feature/your-feature
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_build_disorder_atlas.py

# Run with coverage
pytest --cov=scripts/ tests/

# Run specific test
pytest tests/test_build_disorder_atlas.py::test_load_fasta_sequences
```

### Testing Your Changes

```bash
# Dry run Snakemake to check workflow
snakemake --dry-run --cores 1

# Run specific rule
snakemake --cores 4 build_atlas

# Validate outputs
python3 -c "import pandas as pd; df = pd.read_csv('results/atlas/TP53_atlas_with_clinvar.tsv', sep='\t'); print(f'Loaded {len(df)} residues'); print(df.isnull().sum())"
```

---

## Reporting Security Issues

🔒 **Do not open a public issue for security vulnerabilities.**

Instead, email: s.koganti@[institution].edu with:
- Description of vulnerability
- Steps to reproduce
- Potential impact
- Suggested fix (if any)

---

## Contributor Recognition

All contributors are recognized in:
1. GitHub Contributors page
2. CONTRIBUTORS.md file
3. Acknowledgments section of publications
4. Release notes

---

## Questions?

- 📧 Email: s.koganti@[institution].edu
- 💬 Discussions: GitHub Discussions tab
- 🐛 Issues: GitHub Issues for bugs
- 📚 Docs: See README.md and SUPPLEMENTARY_INFORMATION.md

---

## License

By contributing to this project, you agree that your contributions will be licensed under the same license as the project (CC-BY-4.0 for data, MIT for code).

---

## Acknowledgments

This contributing guide is inspired by best practices from major open-source projects. Thank you to all contributors who help make this resource better!

---

**Last updated:** August 2026  
**Maintainers:** Siva Koganti (@SivaKoganti)
