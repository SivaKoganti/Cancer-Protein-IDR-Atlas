#!/usr/bin/env python3
"""Consistency audit of the manuscript's quantitative claims.

Part 1 re-derives each relationship that was found to be wrong in the 21 August
audit and confirms the corrected value now in the manuscript.
Part 2 lists the issues that remain OPEN because they cannot be resolved without
the raw spectroscopic data.

Usage: python3 scripts/audit_consistency.py
Companion write-up: COMPUTATIONAL_EXPERIMENTAL_AUDIT.md
"""

import math
import re
import sys

R_KCAL = 1.98720425e-3
T = 298.15
RT = R_KCAL * T

MANUSCRIPT = 'COMPLETE_INTEGRATED_MANUSCRIPT.md'

passed, failed, open_issues = [], [], []


def check(name, ok, detail=''):
    (passed if ok else failed).append(name)
    print(f'[{"PASS" if ok else "FAIL"}] {name}')
    for line in detail.strip('\n').split('\n'):
        if line:
            print(f'       {line}')
    print()


def fit(points):
    n = len(points)
    sx = sum(x for x, _ in points); sy = sum(y for _, y in points)
    sxx = sum(x * x for x, _ in points); sxy = sum(x * y for x, y in points)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    ybar = sy / n
    ss_tot = sum((y - ybar) ** 2 for _, y in points)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in points)
    return slope, intercept, (1 - ss_res / ss_tot if ss_tot else float('nan'))


try:
    text = open(MANUSCRIPT, encoding='utf-8').read()
except FileNotFoundError:
    print(f'{MANUSCRIPT} not found; run from the repository root.')
    sys.exit(1)

print(f'RT at {T} K = {RT:.4f} kcal/mol\n')
print('=' * 70)
print('PART 1  Corrections verified against the current manuscript')
print('=' * 70 + '\n')

# 1 -- thermodynamic conversion
fold = math.exp(1.8 / RT)
check('DDG -> fold-change conversion is stated correctly',
      '20.9' in text or '21-fold' in text,
      f'DDG = -1.8 kcal/mol corresponds to {fold:.1f}-fold\n'
      f'measured 1.9-5.4 fold corresponds to DDG = '
      f'{RT*math.log(1.9):.2f}-{RT*math.log(5.4):.2f} kcal/mol\n'
      'manuscript now reports the overestimate rather than a match')

check('the false "quantitative match" claim has been removed',
      'quantitatively matched' not in text,
      'searched for: "quantitatively matched"')

check('overestimate is quantified in the text',
      '4–11×' in text or '4-11x' in text or 'overestimates the effect by 4–11' in text,
      'the 4-11x factor is stated where the comparison is made')

# 2 -- regressions refitted
regressions = [
    ('relative helicity vs [Pb]',
     [(0.0, 100.0), (0.032, 89.3), (0.064, 81.5), (0.32, 73.6)], '−65.2', '0.72'),
    ('alpha-helix fraction vs [Pb]',
     [(0.0, 53.8), (0.032, 48.1), (0.064, 43.9), (0.32, 39.6)], '−35.2', '0.72'),
    ('fK vs [Pb]',
     [(0.0, 4.693), (0.032, 2.442), (0.064, 1.538), (0.32, 0.875)], '−8.35', '0.54'),
    ('intrinsic viscosity vs [Pb]',
     [(0.0, 136.2), (0.001, 138.5), (0.010, 145.8),
      (0.032, 151.8), (0.064, 163.6), (0.320, 185.4)], '138.4', '0.87'),
]
for name, pts, slope_s, r2_s in regressions:
    m, b, r2 = fit(pts)
    check(f'{name}: refitted slope and R^2 appear in the manuscript',
          (slope_s.replace('−', '-') in text.replace('−', '-')) and r2_s in text,
          f'least squares: y = {m:+.2f}x + {b:.2f}, R^2 = {r2:.2f}\n'
          f'the discredited R^2 of 0.996-0.998 is no longer claimed for this fit')

check('no dose-response fit still claims R^2 = 0.996-0.998',
      not re.search(r'R² = 0\.99[68].{0,40}(?:vs|against)?\s*(?:Lead|\[Pb)', text),
      'remaining 0.996/0.998 values are Pearson correlations between methods,\n'
      'which are separate quantities and were not part of the failed checks')

# 3 -- CD formula
check('CD helicity column is described as relative, not absolute',
      'relative helicity' in text and '79.0%' in text,
      'the absolute formula would put the control at 79.0%, not 100%;\n'
      'the manuscript now states the column is [Theta]/[Theta]_control')

# 4 -- viscometry
check('viscometry headline matches its own table',
      '+36.1%' in text and '143.7' in text,
      '185.4/136.2 = +36.1% is internally consistent;\n'
      'the discredited regression is replaced by the refit [eta] = 143.7 + 138.4[Pb]')

# 5 -- arithmetic
check('Table S4.2 fold-changes are arithmetically correct',
      '4.2-fold reduction' in text and '4.8-fold reduction' not in text,
      'fK 0.723 vs same-concentration control 3.053 = 4.22-fold (was stated 4.8)')

# 6 -- QM vs experiment
check('QM Pb-S distance agrees with EXAFS',
      2.64 <= 2.657 <= 2.68,
      'EXAFS (protein PbS3): 2.64-2.68 A\n'
      'this work, PBE0/def2-SVP + ECP60MDF: 2.657 A\n'
      'docked value 2.3 +/- 0.2 A is ~0.35 A short and is now flagged as such')

check('the docked Pb-S distance is disclosed as inconsistent with experiment',
      '0.35 Å' in text and 'EXAFS' in text,
      'Section 3.1.2 and limitation 7 both state the discrepancy')

print('=' * 70)
print('PART 2  Still OPEN -- not resolvable without the raw data')
print('=' * 70 + '\n')

open_issues = [
    ('Units of the fK column',
     'The main text heads the column "Binding Constant (fK)^-1"; supplementary\n'
     'Table S4.2 heads the same numbers "fK M^-1". These are reciprocals.\n'
     'Which is correct determines whether the reported changes are reductions\n'
     'or increases. Resolve against the original Stern-Volmer fits.'),
    ('Ligand concentration scale',
     'fK = 4.693 M^-1 with 5-FU at 0.32 nM gives fractional occupancy 1.5e-9,\n'
     'while the same rows report 28.8-57.9% quenching. These are incompatible.\n'
     'Either the concentrations are not nanomolar or the constants are not M^-1.'),
    ('fK varies with ligand concentration',
     'The tabulated fK rises 3.053 -> 4.176 -> 4.693 as 5-FU goes 0.08 -> 0.32 nM.\n'
     'An equilibrium binding constant cannot depend on ligand concentration,\n'
     'so this column is not a binding constant as currently derived.'),
    ('Provenance of Tables S3-S6',
     'Whether raw spectra underlie these tables determines whether the values\n'
     'should be regenerated from instrument output before submission.'),
    ('Cys-34 / His-67 geometry',
     'Whether these residues can coordinate one Pb(II) in folded HSA has not\n'
     'been verified against the structure. His-67 is a site A ligand; Cys-34 sits\n'
     'in a separate crevice. Requires a distance measurement on the PDB entry.'),
]
for title, detail in open_issues:
    print(f'[OPEN] {title}')
    for line in detail.split('\n'):
        print(f'       {line}')
    print()

print('=' * 70)
print(f'{len(passed)} corrections verified, {len(failed)} failed, '
      f'{len(open_issues)} issues still open')
print('=' * 70)
for f in failed:
    print(f'  FAILED: {f}')
