#!/usr/bin/env python3
"""Recompute every claimed computational/experimental agreement in the manuscript.

Each check takes numbers printed in COMPLETE_INTEGRATED_MANUSCRIPT.md and asks
whether the stated relationship between them holds. Nothing here is fitted or
estimated; it is arithmetic on the manuscript's own values.

Usage: python3 scripts/audit_consistency.py
Companion write-up: COMPUTATIONAL_EXPERIMENTAL_AUDIT.md
"""

import math

R_KCAL = 1.98720425e-3          # kcal / (mol K)
T = 298.15
RT = R_KCAL * T

FAILURES = []


def check(name, ok, detail):
    tag = 'PASS' if ok else 'FAIL'
    if not ok:
        FAILURES.append(name)
    print(f'[{tag}] {name}')
    for line in detail.strip('\n').split('\n'):
        print(f'       {line}')
    print()


def r_squared(points, slope, intercept):
    """R^2 of a stated line against the points it is supposed to describe."""
    ybar = sum(y for _, y in points) / len(points)
    ss_tot = sum((y - ybar) ** 2 for _, y in points)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in points)
    return 1 - ss_res / ss_tot if ss_tot else float('nan')


print(f'RT at {T} K = {RT:.4f} kcal/mol\n')

# ---------------------------------------------------------------- 1
ddg = 1.8
fold_from_ddg = math.exp(ddg / RT)
check(
    'DDG of -1.8 kcal/mol equals a 2-4 fold affinity reduction',
    2 <= fold_from_ddg <= 4,
    f'DDG = -{ddg} kcal/mol  ->  {fold_from_ddg:.1f}-fold\n'
    f'2-fold -> DDG = {RT*math.log(2):.2f} ; 4-fold -> DDG = {RT*math.log(4):.2f} kcal/mol\n'
    f'docking overpredicts by ~{fold_from_ddg/3:.0f}x',
)

# ---------------------------------------------------------------- 2
K, L = 4.693, 0.32e-9
theta = K * L / (1 + K * L)
check(
    'fK = 4.693 /M with 5-FU at 0.32 nM can produce ~50% quenching',
    theta > 0.1,
    f'fractional occupancy = {theta:.2e}\n'
    f'reported quenching at these conditions: 28.8-57.9 %\n'
    f'K needed for 50% occupancy at 0.32 nM: {1/L:.1e} /M '
    f'(HSA-drug constants are 1e3-1e6 /M)',
)

# ---------------------------------------------------------------- 3
fk_vs_ligand = [(0.08, 3.053), (0.16, 4.176), (0.32, 4.693)]
spread = fk_vs_ligand[-1][1] / fk_vs_ligand[0][1]
check(
    'fK is constant with respect to ligand concentration',
    spread < 1.05,
    'fK: ' + ', '.join(f'{c} nM -> {k}' for c, k in fk_vs_ligand) + '\n'
    f'varies by {spread:.2f}x; an equilibrium constant cannot depend on [ligand]',
)

# ---------------------------------------------------------------- 4
check(
    'the fK column carries a consistent unit label',
    False,
    'main text line 304 : "Binding Constant (fK)^-1"\n'
    'supplement S4.2    : "fK  M^-1"\n'
    'identical numbers under reciprocal labels; direction of change depends on which',
)

# ---------------------------------------------------------------- 5..8
regressions = [
    ('Helix % vs [Pb]  y = -6.78x + 88.4', -6.78, 88.4, 0.998,
     [(0.032, 89.3), (0.064, 81.5), (0.32, 73.6)]),
    ('alpha-helix % vs [Pb]  y = -6.78x + 52.1', -6.78, 52.1, 0.996,
     [(0.032, 48.1), (0.064, 43.9), (0.32, 39.6)]),
    ('fK vs [Pb]  y = -13.4x + 4.69', -13.4, 4.69, 0.998,
     [(0.0, 4.693), (0.032, 2.442), (0.064, 1.538), (0.32, 0.875)]),
]
for name, m, b, claimed, pts in regressions:
    actual = r_squared(pts, m, b)
    check(
        f'{name} achieves its claimed R^2 = {claimed}',
        actual >= claimed - 0.05,
        f'actual R^2 against the tabulated points = {actual:.3f}\n'
        + ('negative R^2: fits worse than a horizontal line'
           if actual < 0 else ''),
    )

check(
    'fK regression stays physical across the studied [Pb] range',
    (4.69 - 13.4 * 0.64) > 0,
    f'fit at [Pb] = 0.64 mM -> {4.69 - 13.4*0.64:.2f} (negative binding constant)',
)

# ---------------------------------------------------------------- 9
theta222 = {'Control': -31200, 'Pb 0.032': -27850, 'Pb 0.064': -25420,
            'Pb 0.32': -22950, '5-FU': -30700, 'Pb+5FU': -21480}
stated = {'Control': 100.0, 'Pb 0.032': 89.3, 'Pb 0.064': 81.5,
          'Pb 0.32': 73.6, '5-FU': 98.4, 'Pb+5FU': 68.8}
formula = {k: (v / -39500) * 100 for k, v in theta222.items()}
worst = max(abs(stated[k] - formula[k]) for k in stated)
check(
    'CD helicity matches its stated formula ([Theta]222 / -39500) x 100',
    worst < 2.0,
    '\n'.join(f'{k:<10} stated {stated[k]:>6.1f} %  formula {formula[k]:>6.1f} %'
              for k in theta222) + '\n'
    'the column is [Theta]/[Theta]_control (relative), not helical content',
)

# ---------------------------------------------------------------- 10
eta0, slope_eta = 136.2, 96.8
pct = {c: slope_eta * c / eta0 * 100 for c in (0.032, 0.064, 0.32, 0.64)}
check(
    'viscometry headline +36.1% follows from [eta] = 136.2 + 96.8[Pb]',
    any(abs(v - 36.1) < 1.0 for v in pct.values()),
    '\n'.join(f'[Pb] = {c} mM -> +{v:.1f} %' for c, v in pct.items()) + '\n'
    f'+36.1 % would require [Pb] = {0.361*eta0/slope_eta:.2f} mM, not a studied concentration',
)

# ---------------------------------------------------------------- 11
ref_same = {0.08: 3.053, 0.16: 4.176, 0.32: 4.693}
rows = [(0.08, 0.723, 4.8), (0.16, 1.179, 3.5), (0.32, 2.442, 1.9),
        (0.32, 1.538, 3.1), (0.32, 0.875, 5.4)]
bad = [(v, c) for conc, v, c in rows
       if min(abs(ref_same[conc] / v - c), abs(4.693 / v - c)) > 0.15]
check(
    'Table S4.2 fold-change arithmetic is correct',
    not bad,
    '\n'.join(f'fK = {v}: claimed {c}-fold, computed '
              f'{ref_same[conc]/v:.2f} (same-conc) / {4.693/v:.2f} (vs 4.693)'
              for conc, v, c in rows),
)

# ---------------------------------------------------------------- QM vs experiment
print('=' * 68)
print('QM performed in this session, against the published experimental literature')
print('=' * 68)
print(f'  {"EXAFS, protein/peptide PbS3":<44} 2.64-2.68 A')
print(f'  {"this work, PBE0/def2-SVP + ECP60MDF":<44} 2.657 A     <- agrees')
print(f'  {"manuscript Table 3.1.1 (AutoDock Vina)":<44} 2.3 +/- 0.2 A  <- 0.35 A short')
print(f'  {"Mulliken charge on Pb (first-shell model)":<44} +0.744 e (formal +2)')

print()
print('=' * 68)
print(f'{len(FAILURES)} of 11 checks FAILED')
print('=' * 68)
for f in FAILURES:
    print(f'  - {f}')
