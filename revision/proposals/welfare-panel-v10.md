# Welfare panel --- the 5x3 matrix (generated preview)

Generated from `runs/matrix_5x3-v10-*/manifest.toml` (schema v2) by
`welfare-panel.py`; no number is hand-typed. PREVIEW ONLY: the exact money-metric
is reconstructed from manifest fields, valid because the consumption price index is
the numeraire (`Pi(p) = 1` on every solve path). The kernel-side implementation
supersedes it and is the only citable route. Percent against the reference
continuation; `BETA` is the sectoral cell at `eta_s = 0.5`.

## Table 1 --- the three indices

`(b)` and `(c)` are on the `c+g` bundle with `lambda = 1`, `sigma = 0.9`; they are
reported only where prices are pinned. `(c.C)` is the same leisure correction on
the consumption base, computable for all fifteen cells.

| Cell | dY | dL | dP | (a) Tornqvist | (a) exact EV | (b) c+g | (c) c+g def | (c) c+g leak | (c.C) def | (c.C) leak |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | -0.011 | -0.000 | +0.633 | +0.096 | +0.142 | n/a | n/a | n/a | +0.142 | +0.161 |
| `BF-F2` | -0.016 | -0.000 | +0.718 | -1.980 | -1.980 | n/a | n/a | n/a | -1.980 | -2.250 |
| `BF-F3` | -0.016 | -0.000 | +0.718 | -1.980 | -1.980 | n/a | n/a | n/a | -1.980 | -2.250 |
| `ALPHA-F1` | +0.000 | +0.000 | +0.000 | +0.143 | +0.100 | +0.076 | +0.078 | +0.076 | +0.100 | +0.113 |
| `ALPHA-F2` | +0.000 | +0.000 | +0.000 | -1.823 | -1.823 | +0.078 | -0.101 | +0.078 | -1.823 | -2.071 |
| `ALPHA-F3` | +0.000 | +0.000 | +0.000 | -1.823 | -1.823 | +0.078 | -0.101 | +0.078 | -1.823 | -2.071 |
| `BETA-F1` | +0.122 | +0.119 | +0.270 | +0.423 | +0.420 | n/a | n/a | n/a | +0.268 | +0.304 |
| `BETA-F2` | +0.129 | +0.126 | +0.287 | -1.571 | -1.571 | n/a | n/a | n/a | -1.731 | -1.967 |
| `BETA-F3` | +0.129 | +0.126 | +0.287 | -1.571 | -1.571 | n/a | n/a | n/a | -1.731 | -1.967 |
| `GAMMA-F1` | -0.097 | -0.097 | +0.000 | -0.081 | -0.124 | -0.095 | +0.000 | +0.013 | -0.000 | -0.000 |
| `GAMMA-F2` | +0.126 | +0.126 | +0.000 | -1.533 | -1.533 | +0.299 | +0.000 | +0.160 | -1.694 | -1.924 |
| `GAMMA-F3` | +1.780 | +1.780 | -0.000 | +2.264 | +2.264 | +3.199 | +1.331 | +1.234 | -0.000 | -0.000 |
| `DELTA-F1` | -0.097 | -0.097 | +0.000 | -0.081 | -0.124 | -0.095 | +0.000 | +0.013 | -0.000 | -0.000 |
| `DELTA-F2` | +0.126 | +0.126 | +0.000 | -1.533 | -1.533 | +0.299 | +0.000 | +0.160 | -1.694 | -1.924 |
| `DELTA-F3` | +1.780 | +1.780 | +0.000 | +2.264 | +2.264 | +3.199 | +1.331 | +1.234 | -0.000 | -0.000 |

## Table 2 --- the public-valuation weight (pinned cells)

`(Eh + lam*g1)/(C0 + lam*G0) - 1`. `lam = 0` is the consumption measure of Table 1.

| Cell | lam=0 | lam=0.25 | lam=0.5 | lam=0.75 | lam=1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `ALPHA-F1` | +0.100 | +0.093 | +0.086 | +0.081 | +0.076 |
| `ALPHA-F2` | -1.823 | -1.245 | -0.745 | -0.308 | +0.078 |
| `ALPHA-F3` | -1.823 | -1.245 | -0.745 | -0.308 | +0.078 |
| `GAMMA-F1` | -0.124 | -0.115 | -0.107 | -0.100 | -0.095 |
| `GAMMA-F2` | -1.533 | -0.977 | -0.495 | -0.073 | +0.299 |
| `GAMMA-F3` | +2.264 | +2.548 | +2.794 | +3.009 | +3.199 |
| `DELTA-F1` | -0.124 | -0.115 | -0.107 | -0.100 | -0.095 |
| `DELTA-F2` | -1.533 | -0.977 | -0.495 | -0.073 | +0.299 |
| `DELTA-F3` | +2.264 | +2.548 | +2.794 | +3.009 | +3.199 |

## Table 3 --- the saving convention, consumption base (all fifteen cells)

`D` reads `sE` as deferred consumption, `L` as a leak. The two bases differ by
`E0 - Eh0 = s*E0`, so `L - D = s*(dE - w*dN)/Eh0`, not `s*w*dN/Eh0`.

| Cell | dN | (c.C) deferred | (c.C) leak | difference |
| --- | ---: | ---: | ---: | ---: |
| `BF-F1` | -0.000 | +0.142 | +0.161 | +0.019 |
| `BF-F2` | -0.000 | -1.980 | -2.250 | -0.270 |
| `BF-F3` | -0.000 | -1.980 | -2.250 | -0.270 |
| `ALPHA-F1` | +0.000 | +0.100 | +0.113 | +0.014 |
| `ALPHA-F2` | +0.000 | -1.823 | -2.071 | -0.248 |
| `ALPHA-F3` | +0.000 | -1.823 | -2.071 | -0.248 |
| `BETA-F1` | +0.119 | +0.268 | +0.304 | +0.036 |
| `BETA-F2` | +0.126 | -1.731 | -1.967 | -0.236 |
| `BETA-F3` | +0.126 | -1.731 | -1.967 | -0.236 |
| `GAMMA-F1` | -0.097 | -0.000 | -0.000 | -0.000 |
| `GAMMA-F2` | +0.126 | -1.694 | -1.924 | -0.231 |
| `GAMMA-F3` | +1.780 | -0.000 | -0.000 | -0.000 |
| `DELTA-F1` | -0.097 | -0.000 | -0.000 | -0.000 |
| `DELTA-F2` | +0.126 | -1.694 | -1.924 | -0.231 |
| `DELTA-F3` | +1.780 | -0.000 | -0.000 | -0.000 |

## Blocked and open

- `(b)` and `(c)` on the `c+g` base are `n/a` for the six price-moving cells
  (`BF-F*`, `BETA-F*-etas05`): the bundle deflator is not in the manifests.
- `(c.C)` for `BF` values leisure at the wage-bill-weighted aggregate wage; it
  should use the binding sectoral wage.
- All leisure-augmented columns are first-order. The second-order term needs the
  curvature specification (D5).
- `sigma` is irrelevant wherever prices are pinned: `(b)` is `sigma`-free in 9 of 15 cells.
