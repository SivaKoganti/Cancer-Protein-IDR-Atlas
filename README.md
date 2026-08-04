# Cancer Protein IDR & LLPS Atlas

This repository contains the code and pipeline to generate intrinsic disorder region (IDR), low-complexity region (LCR), and liquid–liquid phase separation (LLPS) propensity annotations for cancer-associated human proteins and to integrate ClinVar variants and evolutionary analyses.

This project is an MVP scaffold and includes:

- Snakemake pipeline to download data and run predictors (IUPred2A, SEG, PLAAC) and map ClinVar variants.
- Dockerfile and requirements for reproducible execution.
- Scripts to run batch predictors and aggregate results.
- A curated starter list of 50 high-confidence Cancer Gene Census genes for the MVP.

Next steps (automated after this commit):
1. Run `snakemake` (or `docker build` + `snakemake`) to download ClinVar and UniProt proteome and extract the selected proteins.
2. Install predictor binaries (IUPred2A, PLAAC, SEG) in the container or system path.
3. Run the pipeline to create per-protein annotation tables and sample plots.

If you want me to run the pilot and push results, grant me details on where to store large outputs (Git LFS, S3) or I can upload small sample outputs to this repo.

License: MIT (see LICENSE)
