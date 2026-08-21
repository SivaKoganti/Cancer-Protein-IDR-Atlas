# The Lead–HSA–5-FU manuscript has moved out of this repository

That work is unrelated to the Cancer Protein IDR Atlas and no longer lives here.
It is now a self-contained git repository with its own history, delivered as a
bundle: `lead-hsa-5fu-manuscript.bundle`.

## Extract it

```bash
git clone lead-hsa-5fu-manuscript.bundle lead-hsa-5fu-manuscript
cd lead-hsa-5fu-manuscript
git remote remove origin        # detach from the bundle file
```

You now have an ordinary git repository. To give it a home on GitHub, create an
empty repository there and push:

```bash
git remote add origin git@github.com:<you>/lead-hsa-5fu-manuscript.git
git push -u origin main
```

## What it contains

The manuscript (Markdown source and Word build with figures embedded), three cover
letters, 13 figures at 300 dpi, supplementary Tables S3–S6, the docking inputs and
outputs, and the analysis code. `scripts/audit_consistency.py` recomputes every
quantitative claim the manuscript makes and currently reports 20 checks passed,
0 failed. Its README lists the three items still awaiting the authors.

## Why a bundle

A git bundle is a single file holding a complete repository — history included.
It was used here because this environment cannot create a new GitHub repository,
and the manuscript needed a durable home that shares no history with this project.
Once extracted and pushed, this file can be deleted.
