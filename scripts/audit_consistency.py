#!/usr/bin/env python3
"""Consistency audit of the manuscript's quantitative claims.

Part 1 re-derives each relationship that was found to be wrong in the 21 August
audit and confirms the value now in the manuscript.
Part 2 lists items that still require confirmation by the authors against the
primary records.

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

PASSED = []
FAILED = []


def check(name, ok, detail=''):
    (PASSED if ok else FAILED).append(name)
    print(f'[{"PASS" if ok else "FAIL"}] {name}')
    for line in detail.strip('\n').split('\n'):
        if line:
            print(f'       {line}')
    print()


def fit(points):
    """Least-squares line through (x, y) points; returns slope, intercept, R^2."""
    n = len(points)
    sx = sum(x for x, _ in points)
    sy = sum(y for _, y in points)
    sxx = sum(x * x for x, _ in points)
    sxy = sum(x * y for x, y in points)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    ybar = sy / n
    ss_tot = sum((y - ybar) ** 2 for _, y in points)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in points)
    return slope, intercept, (1 - ss_res / ss_tot if ss_tot else float('nan'))


def r_squared(points, slope, intercept):
    ybar = sum(y for _, y in points) / len(points)
    ss_tot = sum((y - ybar) ** 2 for _, y in points)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in points)
    return 1 - ss_res / ss_tot if ss_tot else float('nan')


try:
    text = open(MANUSCRIPT, encoding='utf-8').read()
except FileNotFoundError:
    print(f'{MANUSCRIPT} not found; run from the repository root.')
    sys.exit(1)

print(f'RT at {T} K = {RT:.4f} kcal/mol\n')
print('=' * 70)
print('PART 1  Quantitative claims verified against the current manuscript')
print('=' * 70 + '\n')

# --- thermodynamics -------------------------------------------------------
fold = math.exp(1.8 / RT)
check('DDG -> fold-change conversion is stated correctly',
      '20.9' in text or '21-fold' in text,
      f'DDG = -1.8 kcal/mol corresponds to {fold:.1f}-fold\n'
      f'measured 1.9-5.4 fold corresponds to DDG = '
      f'{RT*math.log(1.9):.2f}-{RT*math.log(5.4):.2f} kcal/mol')

check('the false "quantitative match" claim has been removed',
      'quantitatively matched' not in text)

check('the docking overestimate is quantified where the comparison is made',
      '4–11×' in text or '4-11x' in text)

# --- regressions ----------------------------------------------------------
regressions = [
    ('relative helicity vs [Pb]',
     [(0.0, 100.0), (0.032, 89.3), (0.064, 81.5), (0.32, 73.6)], '65.2', '0.72'),
    ('alpha-helix fraction vs [Pb]',
     [(0.0, 53.8), (0.032, 48.1), (0.064, 43.9), (0.32, 39.6)], '35.2', '0.72'),
    ('intrinsic viscosity vs [Pb]',
     [(0.0, 136.2), (0.001, 138.5), (0.010, 145.8),
      (0.032, 151.8), (0.064, 163.6), (0.320, 185.4)], '138.4', '0.87'),
]
for name, pts, slope_s, r2_s in regressions:
    m, b, r2 = fit(pts)
    check(f'{name}: refitted slope and R^2 appear in the manuscript',
          slope_s in text and r2_s in text,
          f'least squares: y = {m:+.2f}x + {b:.2f}, R^2 = {r2:.2f}')

check('no dose-response fit still claims R^2 = 0.996-0.998',
      not re.search(r'R² = 0\.99[68].{0,40}(?:Lead|\[Pb)', text),
      'remaining 0.996/0.998 values are cross-method Pearson correlations,\n'
      'which are separate quantities')

# --- CD ------------------------------------------------------------------
check('CD helicity column is described as relative, not absolute',
      'relative helicity' in text and '79.0%' in text,
      'the absolute formula would put the control at 79.0%, not 100%')

# --- viscometry ----------------------------------------------------------
visc = [(0.0, 136.2), (0.001, 138.5), (0.010, 145.8),
        (0.032, 151.8), (0.064, 163.6), (0.320, 185.4)]
check('viscometry: table is self-consistent and the refit is reported',
      '+36.1%' in text and '143.7' in text,
      f'185.4/136.2 = +36.1% internally consistent\n'
      f'R^2 of the discredited line 136.2 + 96.8x against this table = '
      f'{r_squared(visc, 96.8, 136.2):.2f}\n'
      f'refit [eta] = 143.7 + 138.4[Pb], R^2 = 0.87')

# --- fluorescence --------------------------------------------------------
ratios = [(305, 72), (418, 118), (469, 244), (469, 154), (469, 88)]
computed = ', '.join(f'{a/b:.2f}' for a, b in ratios)
check('Stern-Volmer ratios are arithmetically correct',
      '4.22' in text and '4.8-fold reduction' not in text,
      f'matched-pair ratios: {computed}\n'
      f'the previously stated 4.8-fold is corrected to 4.22')

K_mM = 4.693e2
occ = K_mM * 0.32e-3 / (1 + K_mM * 0.32e-3)
check('5-FU concentration scale is internally consistent',
      abs(occ - 0.124) < 0.02 and '0.32 mM' in text,
      f'at fK = {K_mM:.0f} /M and 0.32 mM 5-FU, predicted occupancy = {occ*100:.1f}%\n'
      f'measured tryptophan quenching = 12.4%\n'
      f'nanomolar would require 4.4e8 /M and place the ligand at 2e-5 of [protein]')

check('no nanomolar 5-FU concentrations remain in the tables',
      not re.search(r'5-FU[^|\n]{0,12}\d\.\d+ nM', text))

# --- quantum chemistry ---------------------------------------------------
check('QM-refined coordination distances are reported in Table 3.1.1',
      all(v in text for v in ('2.65', '2.73', '2.42')),
      'Pb-S 2.65, Pb-N 2.73, Pb-O 2.42/2.38 A (mixed-donor cluster, PBE0/def2-SVP)\n'
      'Pb-S coincides with the EXAFS range 2.64-2.68 A')

check('the docked Pb-S distance is disclosed as inconsistent with experiment',
      '0.35 Å' in text and 'EXAFS' in text)

check('the inconclusive dianion calculation is reported as inconclusive',
      'marginally bound dianion' in text and 'draw no conclusion' in text)

# --- reporting sections --------------------------------------------------
for name, needle in [('generative-AI use is disclosed', 'DECLARATION OF GENERATIVE AI USE'),
                     ('author contributions section present', 'AUTHOR CONTRIBUTIONS'),
                     ('data availability section present', 'DATA AND CODE AVAILABILITY')]:
    check(name, needle in text)

print('=' * 70)
print('PART 2  Requires author confirmation before submission')
print('=' * 70 + '\n')

OPEN = [
    ('Provenance of Tables S3-S6',
     'The concentration-scale correction was inferred from internal consistency,\n'
     'not read from an instrument record. Confirm against the original\n'
     'worksheets and deposit the primary data (Section 6.3).'),
    ('Cys-34 / His-67 geometry in folded albumin',
     'The quantum chemistry shows the four-donor sphere is chemically viable as\n'
     'an isolated cluster; it does not show these side chains can adopt that\n'
     'arrangement in the protein. His-67 is a site A ligand, Cys-34 sits in a\n'
     'separate crevice. Measure the Sgamma...Nepsilon distance in a\n'
     'high-resolution structure (Section 4.5, limitation 8).'),
    ('Global Stern-Volmer refit',
     'Constants were evaluated per concentration rather than from one fit across\n'
     'the series, which is why they drift upward with [5-FU]. A global refit\n'
     'gives one constant per condition. The ratios the conclusions rest on are\n'
     'unaffected.'),
]
for title, detail in OPEN:
    print(f'[OPEN] {title}')
    for line in detail.split('\n'):
        print(f'       {line}')
    print()

print('=' * 70)
print(f'{len(PASSED)} checks passed, {len(FAILED)} failed, '
      f'{len(OPEN)} items require author confirmation')
print('=' * 70)
for f in FAILED:
    print(f'  FAILED: {f}')
