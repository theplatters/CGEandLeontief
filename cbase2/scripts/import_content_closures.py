#!/usr/bin/env python3
"""Import content of the programme across the closure matrix --- NOT INCLUDED in the
manuscript. The comparison table this script writes was withdrawn on JK's instruction
(the full-model multiplier is not a measure of the programme's import content: its
numerator grows with the response, so it ranks closures by output volume rather than by
import intensity). Appendix C remarks why; the script and its table are kept as the
measurement behind that remark.

Writes revision/tables/table_import-content-closures.tex.

The programme's own spending is a cost-share object and is closure-invariant;
what the closures change is how much output (and therefore how many imported
intermediates) the programme actually sets in motion. This script puts the two
side by side:

  benchmark   the programme's domestic bill through the domestic Leontief
              inverse (no induced household demand) --- the figures of
              Appendix C.9 and of composition.py
  closures    the net gross-output change of the reported runs, priced at the
              table's sectoral import intensities --- i.e. including the
              induced household demand and the crowding-out of the financing

Sources (all by label or by committed run id, never by position):
  data/I-O_DE2019_formatiert.csv       the 2019 use table
  data/impulses.csv                    the programme (2024 row)
  cbase2/data_processed/AC_domar_weights.csv   baseline gross output, GDP units
  runs/matrix_5x3-v8-<cell>/solution.csv       sectoral solutions (v8 ships them)
  runs/matrix_5x3-v11-<cell>/manifest.toml     the reported headline metrics

Canaries (the build aborts on either):
  1. v8 and v11 agree on gdp_rel to 1e-9 for every cell --- so the family that
     ships solutions is the family the paper reports.
  2. sum_i (v_i/y_i) dq_i reproduces gdp_rel to 1e-8 --- the GDP identity, which
     fixes the units of the solution vectors.

The comparison is gross versus net by construction, and the caption says so:
the benchmark is the gross output the programme's own domestic spending calls
forth, the cells are the net change the model produces once the financing and
the induced household demand are included. Import content is the imports
generated, as a share of the programme; the bill's own import share (17.8 per
cent) is common to every row and is not repeated.

Usage:  python3 cbase2/scripts/import_content_closures.py [--check]
"""
from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
TABLE = DATA / "I-O_DE2019_formatiert.csv"
IMPULSES = DATA / "impulses.csv"
LAMBDA = ROOT / "cbase2/data_processed/AC_domar_weights.csv"
RUNS = ROOT / "runs"
OUT = ROOT / "revision/tables/table_import-content-closures.tex"

G0 = 40_300.0          # programme, mn EUR, 2019 prices
GDP = 3_027_818.0      # GDP at basic prices, mn EUR (the model's numeraire = 1)
G = G0 / GDP           # the programme in GDP units

CELLS = [f"{c}-{f}" for c in ("BF", "ALPHA", "BETA", "GAMMA", "DELTA")
         for f in ("F1", "F2", "F3")]

# grouped presentation: the equivalence classes the runs actually show
GROUPS = [
    ("BF, F1", ["BF-F1"]),
    ("BF, F2--F3", ["BF-F2", "BF-F3"]),
    ("ALPHA and BETA, F1", ["ALPHA-F1", "BETA-F1"]),
    ("ALPHA and BETA, F2--F3", ["ALPHA-F2", "ALPHA-F3", "BETA-F2", "BETA-F3"]),
    ("GAMMA and DELTA, F1", ["GAMMA-F1", "DELTA-F1"]),
    ("GAMMA and DELTA, F2", ["GAMMA-F2", "DELTA-F2"]),
    ("GAMMA and DELTA, F3", ["GAMMA-F3", "DELTA-F3"]),
]


def _num(cell: str) -> float:
    s = cell.strip()
    return 0.0 if s in {"-", "", "x"} else float(s)


def read_table():
    rows = list(csv.reader(open(TABLE, encoding="utf-8-sig"), delimiter=";", quotechar='"'))
    header, body = rows[0][1:], rows[1:]
    labels = [r[0].strip() for r in body]
    M = np.array([[_num(c) for c in r[1:]] for r in body], dtype=float)

    def row(prefix: str) -> int:
        hits = [k for k, l in enumerate(labels) if l.startswith(prefix)]
        if not hits:
            raise KeyError(prefix)
        return hits[0]

    def col(prefix: str) -> int:
        hits = [k for k, c in enumerate(header) if c.startswith(prefix)]
        if not hits:
            raise KeyError(prefix)
        return hits[0]

    sectors = list(range(71))
    inv = [col("Anlageinvestitionen f.Ausr\u00fcstungen"), col("Anlageinvestitionen f\u00fcr Bauten")]
    return dict(
        y=M[row("Produktionswert"), sectors],
        Z=M[np.ix_(sectors, sectors)],
        m=M[row("Verwendung der Importe"), sectors],
        v=M[row("Bruttowertsch\u00f6pfung"), sectors],
        inv=inv,
        M=M, row=row, col=col,
    )


def programme() -> np.ndarray:
    imp = pd.read_csv(IMPULSES)
    cols = [c for c in imp.columns if c not in {"", "year", "Unnamed: 0"}][:71]
    row = imp[imp.year == 2024].iloc[0]
    psi = np.array([float(row[c]) for c in cols])
    return G0 * psi / psi.sum()


def metric(run: str, key: str = "gdp_rel") -> float:
    p = RUNS / run / "manifest.toml"
    return float(tomllib.load(open(p, "rb"))["metrics"][key])


def solution(run: str) -> np.ndarray:
    return pd.read_csv(RUNS / run / "solution.csv").quantity.to_numpy()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="abort on a failed canary")
    args = ap.parse_args()

    d = read_table()
    y, Z, m, v = d["y"], d["Z"], d["m"], d["v"]
    P = Z / np.where(y == 0, 1.0, y)
    g = programme()

    # --- bill at purchasers' prices (same convention as composition.py)
    dom_inv = d["M"][d["row"]("Gesamte Verwendung der inl\u00e4ndisch"), d["inv"]].sum()
    imp_inv = d["M"][d["row"]("Verwendung der Importe"), d["inv"]].sum()
    tax_inv = d["M"][d["row"]("G\u00fctersteuern"), d["inv"]].sum()
    tot_inv = dom_inv + imp_inv + tax_inv
    sh_imp = imp_inv / tot_inv

    # --- benchmark: the programme's domestic bill through the Leontief inverse
    g_dom = (dom_inv / tot_inv) * g
    x = np.linalg.solve(np.eye(71) - P, g_dom)
    imp_bench = float((np.where(y > 0, m / y, 0.0) * x).sum())
    bench = dict(bill=100 * sh_imp, interm=100 * imp_bench / G0,
                 total=100 * (sh_imp * G0 + imp_bench) / G0, gross=x.sum() / G0)

    # --- the closures
    lam = pd.read_csv(LAMBDA)["lambda"].astype(float).to_numpy()
    cells = {}
    worst_metric, worst_identity, price_movers = 0.0, 0.0, []
    for c in CELLS:
        dq = solution(f"matrix_5x3-v8-{c}") - lam
        gdp_model = metric(f"matrix_5x3-v11-{c}")
        worst_metric = max(worst_metric, abs(metric(f"matrix_5x3-v8-{c}") - gdp_model))
        gap = abs(float((np.where(y > 0, v / y, 0.0) * dq).sum()) - gdp_model)
        # the GDP identity holds at baseline prices; where the cell moves prices
        # (BF: max |dp| 0.22--0.28) real GDP and baseline-price value added diverge
        if metric(f"matrix_5x3-v11-{c}", "max_abs_price_dev") < 1e-6:
            worst_identity = max(worst_identity, gap)
        else:
            price_movers.append((c, gap))
        interm = float((np.where(y > 0, m / y, 0.0) * dq).sum()) / G
        cells[c] = dict(interm=100 * interm, total=100 * (sh_imp + interm),
                        gross=float(dq.sum()) / G)

    print(f"benchmark   bill {bench['bill']:5.2f} % | intermediates {bench['interm']:5.2f} %"
          f" | total {bench['total']:5.2f} % | output {bench['gross']:5.2f} xG0")
    for c in CELLS:
        r = cells[c]
        print(f"{c:12s} intermediates {r['interm']:5.2f} % | total {r['total']:5.2f} %"
              f" | output {r['gross']:5.2f} xG0")
    print(f"canary 1  max |gdp_rel(v8) - gdp_rel(v11)| = {worst_metric:.2e}")
    print(f"canary 2  max |sum (v/y) dq - gdp_rel|      = {worst_identity:.2e}  (price-stable cells)")
    for c, gap in price_movers:
        print(f"          {c}: identity gap {gap:.2e} -- prices move in this cell, so real GDP and"
              f" baseline-price value added diverge (expected)")

    if args.check and (worst_metric > 1e-9 or worst_identity > 1e-8):
        print("FAILED: a canary is out of tolerance", file=sys.stderr)
        return 1

    # --- grouped rows: a single value where the class coincides, else a range
    def fmt(group: str, key: str, nd: int = 1, tol: float = 0.05) -> str:
        vals = [cells[c][key] for c in group]
        lo, hi = min(vals), max(vals)
        return f"{lo:.{nd}f}" if abs(hi - lo) < tol else f"{lo:.{nd}f}--{hi:.{nd}f}"

    rows = [
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{lrrr}",
        r"\toprule",
        r" & imports & total import & induced \\",
        r" & in production & content & gross output \\",
        r" & (\% of $G_0$) & (\% of $G_0$) & ($\times G_0$) \\",
        r"\midrule",
        r"\multicolumn{4}{l}{\emph{A. The programme's own spending: cost shares, no induced demand}}\\",
        rf"programme's domestic bill & ${bench['interm']:.1f}$ & ${bench['total']:.1f}$ & ${bench['gross']:.2f}$ \\",
        r"\midrule",
        r"\multicolumn{4}{l}{\emph{B. Including the model's induced demand and the financing's crowding-out}}\\",
    ]
    for label, group in GROUPS:
        rows.append(rf"{label} & ${fmt(group, 'interm')}$ & ${fmt(group, 'total')}$ & "
                    rf"${fmt(group, 'gross', nd=2, tol=0.005)}$ \\")
    rows += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\caption{The programme's import content, as a cost-share object and as a model outcome.",
        r"Panel~A prices the programme's domestic bill through the domestic Leontief inverse and",
        r"is invariant across the closure matrix; Panel~B reads the same measure off the net",
        r"gross-output change of each reported run, so it includes induced household demand and the",
        r"crowding-out the financing implies. The bill's own import share, $17.8$ per cent of the",
        r"programme, is common to every row and is not repeated; the intermediate-import column adds",
        r"the imports the production response draws in. The two panels are therefore gross and net",
        r"objects and their difference measures how much of the programme's import content the",
        r"financing and the labour-market closure absorb.}",
        r"\label{tab:import-content-closures}",
        r"\end{table}",
    ]

    header = (
        "% GENERATED by cbase2/scripts/import_content_closures.py -- do not edit by hand.\n"
        "% NOT INCLUDED in the manuscript: withdrawn because the full-model multiplier is not a\n"
        "% measure of the programme's import content (see the remark at the end of Appendix C.9).\n"
        "% Source: data/I-O_DE2019_formatiert.csv (label lookup) + data/impulses.csv (2024 row)\n"
        "%   + runs/matrix_5x3-v8-<cell>/solution.csv + runs/matrix_5x3-v11-<cell>/manifest.toml\n"
        f"% Canary 1 on regeneration: max |gdp_rel(v8) - gdp_rel(v11)| = {worst_metric:.2e}\n"
        f"% Canary 2 on regeneration: max |sum (v/y) dq - gdp_rel|      = {worst_identity:.2e}\n"
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(header + "\n".join(rows) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
