#!/usr/bin/env python3
"""Welfare panel for the 5x3 matrix --- regenerate from the v10 run manifests.

    python3 welfare-panel.py            # writes welfare-panel-v10.md next to this file

No number is hand-typed: every value is a manifest field or an arithmetic
combination of manifest fields, with the mapping stated below. This is a
PREVIEW generator: it reconstructs the exact money-metric from manifest fields,
which is possible because the consumption price index is the numeraire
(Pi(p) = 1 on every solve path). The kernel-side implementation (Phase 1 of the
workplan) supersedes it and is the only citable route.

Conventions (decisions D1-D4, 2026-10-16)
-----------------------------------------
E    household resources, balanced public budget (taxes = purchases):
         F1:  E = nominal_gdp - gdp_g + F            (additive demand = 0)
         F2:  E = nominal_gdp - gdp_g + F            (T_G + T(p) = gdp_g)
         F3:  E = nominal_gdp - gdp_g + B_gov + F    (T(p) = 0)
     Eh  = (1-s) E, the consumption block; s from the calibration.
     g1  = gdp_g = dot(p, gov + programme), the public absorption at cell prices.
Bundles (D1: c+g, lambda = 1):
     private+public, consumption base   Vf = Eh + g1,  Vf0 = 0.905776
     private+public, full-income base   Ve = E  + g1,  Ve0 = 1.000000
Measures:
     (a)  Tornqvist quantity index of consumption        = metrics.consumption_rel
     (a') exact CES money-metric of consumption          = dE/E_0
     (b)  final-demand index over c+g                    = Vf/Vf0 - 1
     (c)  leisure-augmented money-metric on the c+g base, first-order:
              D (saving = deferred consumption)         = (dVe - w dN)/Ve0
              L (saving sE = a leak)                    = (dVf - w dN)/Vf0
     (cC) the same on the consumption base (computable for all 15 cells):
              D = (dE - w dN)/E0,   L = (dE - w dN)/Eh0
     w  = metrics.wage (wage-bill-weighted aggregate), dN = employment - 1.
(b) and (c) on the c+g base are reported ONLY for cells whose prices are pinned
(max_abs_price_dev = 0), where every quantity index reduces to its value ratio.
The six price-moving cells need the bundle-specific deflator Pi_{c+g}(p), which
is not in the manifests.
"""
import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(os.path.dirname(__file__), "welfare-panel-v10.md")

TG_BASE = 0.214101
E0 = 1.0 - TG_BASE          # 0.785899
EH0 = 0.691675
C0, G0 = 0.691675, 0.214101
S = 1.0 - EH0 / E0          # calibrated saving rate
LAB = ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
FIN = ["F1", "F2", "F3"]


def cid(lab, fin):
    return f"matrix_5x3-v10-{lab}-{fin}" + ("-etas05" if lab == "BETA" else "")


def row(lab, fin):
    d = tomllib.load(open(f"{REPO}/runs/{cid(lab, fin)}/manifest.toml", "rb"))
    m, g = d["metrics"], d["diagnostics"]
    F, Bg = m["external_transfer"], m.get("programme_financing", 0.0)
    E = m["nominal_gdp"] - g["gdp_g"] + F + (Bg if fin == "F3" else 0.0)
    Eh, g1 = (1 - S) * E, g["gdp_g"]
    w, dN = m["wage"], m["employment"] - 1.0
    Vf, Vf0 = Eh + g1, C0 + G0
    Ve, Ve0 = E + g1, E0 + G0
    return dict(
        id=f"{lab}-{fin}", pin=m["max_abs_price_dev"] < 1e-12,
        dY=100 * m["gdp_rel"], dN=100 * dN, dP=100 * (m["gdp_deflator"] - 1.0),
        tor=100 * m["consumption_rel"], ev=100 * (E / E0 - 1.0),
        fd=100 * (Vf / Vf0 - 1.0), cD=100 * ((Ve - Ve0) - w * dN) / Ve0,
        cL=100 * ((Vf - Vf0) - w * dN) / Vf0,
        cCD=100 * ((E - E0) - w * dN) / E0, cCL=100 * ((E - E0) - w * dN) / EH0,
        lam=lambda l, Eh=Eh, g1=g1: 100 * ((Eh + l * g1) / (C0 + l * G0) - 1.0),
    )


def f(x):
    return "n/a" if x is None else f"{x:+.3f}"


def build():
    R = [row(l, k) for l in LAB for k in FIN]
    L = []
    A = L.append
    A("# Welfare panel --- the 5x3 matrix (generated preview)")
    A("")
    A("Generated from `runs/matrix_5x3-v10-*/manifest.toml` (schema v2) by")
    A("`welfare-panel.py`; no number is hand-typed. PREVIEW ONLY: the exact money-metric")
    A("is reconstructed from manifest fields, valid because the consumption price index is")
    A("the numeraire (`Pi(p) = 1` on every solve path). The kernel-side implementation")
    A("supersedes it and is the only citable route. Percent against the reference")
    A("continuation; `BETA` is the sectoral cell at `eta_s = 0.5`.")
    A("")
    A("## Table 1 --- the three indices")
    A("")
    A("`(b)` and `(c)` are on the `c+g` bundle with `lambda = 1`, `sigma = 0.9`; they are")
    A("reported only where prices are pinned. `(c.C)` is the same leisure correction on")
    A("the consumption base, computable for all fifteen cells.")
    A("")
    A("| Cell | dY | dL | dP | (a) Tornqvist | (a) exact EV | (b) c+g | (c) c+g def | (c) c+g leak | (c.C) def | (c.C) leak |")
    A("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in R:
        A(f"| `{r['id']}` | {r['dY']:+.3f} | {r['dN']:+.3f} | {r['dP']:+.3f} | {r['tor']:+.3f} | "
          f"{r['ev']:+.3f} | {f(r['fd'] if r['pin'] else None)} | "
          f"{f(r['cD'] if r['pin'] else None)} | {f(r['cL'] if r['pin'] else None)} | "
          f"{r['cCD']:+.3f} | {r['cCL']:+.3f} |")
    A("")
    A("## Table 2 --- the public-valuation weight (pinned cells)")
    A("")
    A("`(Eh + lam*g1)/(C0 + lam*G0) - 1`. `lam = 0` is the consumption measure of Table 1.")
    A("")
    A("| Cell | lam=0 | lam=0.25 | lam=0.5 | lam=0.75 | lam=1 |")
    A("| --- | ---: | ---: | ---: | ---: | ---: |")
    for r in R:
        if r["pin"]:
            A(f"| `{r['id']}` | " + " | ".join(f"{r['lam'](l):+.3f}" for l in (0, .25, .5, .75, 1))
              + " |")
    A("")
    A("## Table 3 --- the saving convention, consumption base (all fifteen cells)")
    A("")
    A("`D` reads `sE` as deferred consumption, `L` as a leak. The two bases differ by")
    A("`E0 - Eh0 = s*E0`, so `L - D = s*(dE - w*dN)/Eh0`, not `s*w*dN/Eh0`.")
    A("")
    A("| Cell | dN | (c.C) deferred | (c.C) leak | difference |")
    A("| --- | ---: | ---: | ---: | ---: |")
    for r in R:
        A(f"| `{r['id']}` | {r['dN']:+.3f} | {r['cCD']:+.3f} | {r['cCL']:+.3f} | "
          f"{r['cCL'] - r['cCD']:+.3f} |")
    A("")
    A("## Blocked and open")
    A("")
    A("- `(b)` and `(c)` on the `c+g` base are `n/a` for the six price-moving cells")
    A("  (`BF-F*`, `BETA-F*-etas05`): the bundle deflator is not in the manifests.")
    A("- `(c.C)` for `BF` values leisure at the wage-bill-weighted aggregate wage; it")
    A("  should use the binding sectoral wage.")
    A("- All leisure-augmented columns are first-order. The second-order term needs the")
    A("  curvature specification (D5).")
    A("- `sigma` is irrelevant wherever prices are pinned: `(b)` is `sigma`-free in 9 of 15 cells.")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    open(OUT, "w").write(build())
    print("wrote", OUT)
