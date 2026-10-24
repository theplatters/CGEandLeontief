#!/usr/bin/env python3
"""Composition of the programme's spending --- generator for Appendix C.9.

Produces revision/tables/table_composition.tex, the three-panel table of
Appendix C.9 (label tab:pipeline-composition), and prints the numbers.

Every object is read from the published table by LABEL, never by position:
  data/I-O_DE2019_formatiert.csv   German 2019 use table (semicolon, plain
                                   integers, '-' and 'x' mean zero)
  data/impulses.csv                the source study's investment impulses;
                                   psi = the year-1 (2024) row over the 71
                                   sector columns, renormalised

Conventions (see Appendix C.9 and the model's own code):
  P[u,i] = Z[u,i] / y_i        the input from u per unit of output of i
  x      = (I - P)^{-1} g_dom  the induced gross output (SUPPLY orientation:
                               the model inverts Omega' diag(1-alpha), which is
                               P. Inverting the transposed system inflates x by
                               about a third -- do not.)
Canary: (I - P) y must equal the table's domestic final demand (y - row sums of
Z); it holds to ~1e-12 and is the check that fixes the orientation.

Usage:  python3 cbase2/scripts/composition.py [--check]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
TABLE = DATA / "I-O_DE2019_formatiert.csv"
IMPULSES = DATA / "impulses.csv"
OUT = ROOT / "revision/tables/table_composition.tex"

G0 = 40_300.0  # the programme's size, mn EUR, 2019 prices (the paper's value)


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
    final_use = [
        col("Konsumausgaben der privaten Haushalte"),
        col("Konsumausgaben der privaten Organisationen"),
        col("Konsumausgaben des Staates"),
        col("Anlageinvestitionen f.Ausr\u00fcstungen"),
        col("Anlageinvestitionen f\u00fcr Bauten"),
        col("Vorratsver\u00e4nderungen"),
        col("Exporte"),  # the total export column; 81-84 are its components
    ]
    investment = final_use[3:5]
    return dict(
        y=M[row("Produktionswert"), sectors],
        Z=M[np.ix_(sectors, sectors)],
        m=M[row("Verwendung der Importe"), sectors],
        t=M[row("G\u00fctersteuern"), sectors],
        v=M[row("Bruttowertsch\u00f6pfung"), sectors],
        dom_final=M[row("Gesamte Verwendung der inl\u00e4ndisch"), final_use],
        imp_final=M[row("Verwendung der Importe"), final_use],
        tax_final=M[row("G\u00fctersteuern"), final_use],
        dom_inv=M[row("Gesamte Verwendung der inl\u00e4ndisch"), investment].sum(),
        imp_inv=M[row("Verwendung der Importe"), investment].sum(),
        tax_inv=M[row("G\u00fctersteuern"), investment].sum(),
    )


def programme() -> np.ndarray:
    import pandas as pd

    imp = pd.read_csv(IMPULSES)
    cols = [c for c in imp.columns if c not in {"", "year", "Unnamed: 0"}][:71]
    row = imp[imp.year == 2024].iloc[0]
    psi = np.array([float(row[c]) for c in cols])
    psi = psi / psi.sum()
    return G0 * psi


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="assert the orientation canary")
    args = ap.parse_args()

    d = read_table()
    y, Z, m, t, v = d["y"], d["Z"], d["m"], d["t"], d["v"]
    P = Z / np.where(y == 0, 1.0, y)          # P[u,i] = input u per unit of output of i
    g = programme()

    # --- canary: the orientation is fixed by the kernel's own calibration identity
    resid = (np.eye(71) - P) @ y
    dom_final_row = y - Z.sum(axis=1)
    canary = float(np.abs(resid - dom_final_row).max())
    print(f"canary  max |(I - P) y - (y - Zrow)| = {canary:.3e}")
    if args.check and canary > 1e-8:
        print("FAILED: wrong orientation", file=sys.stderr)
        return 1

    # --- the bill at purchasers' prices
    tot_inv = d["dom_inv"] + d["imp_inv"] + d["tax_inv"]
    sh_imp, sh_tax = d["imp_inv"] / tot_inv, d["tax_inv"] / tot_inv
    bill = dict(imp=G0 * sh_imp, tax=G0 * sh_tax, dom=G0 * (1 - sh_imp - sh_tax))

    # --- the induced output, by cost structure (content measures)
    g_dom = bill["dom"] / G0 * g
    x = np.linalg.solve(np.eye(71) - P, g_dom)
    inten = lambda vec: float((np.where(y > 0, vec / y, 0.0) * x).sum())
    induced = dict(
        dom_int=inten(Z.sum(axis=0)), imp_int=inten(m), tax_int=inten(t), va=inten(v)
    )
    totals = dict(imp=bill["imp"] + induced["imp_int"], tax=bill["tax"] + induced["tax_int"])

    # --- the band: the direct leg at the programme-weighted import intensity
    psi = g / G0
    w_imp = float((psi * np.where(y > 0, m / y, 0.0)).sum())
    total_high = (bill["imp"] + induced["imp_int"]) / G0
    total_low = w_imp + induced["imp_int"] / G0

    print(f"programme            {G0:12,.0f} mn EUR (2019)")
    print(f"bill                 imported {100*sh_imp:5.2f} % | taxes {100*sh_tax:5.2f} % | domestic {100*(1-sh_imp-sh_tax):5.2f} %")
    print(f"induced output       {x.sum():12,.0f} mn = {x.sum()/G0:.3f} x programme")
    for k, val in induced.items():
        print(f"   {k:18s} {val:10,.0f} mn = {100*val/G0:5.2f} % of programme | {100*val/x.sum():5.2f} % of x")
    print(f"import content       {100*total_high:5.2f} % of the programme (bottom-up {100*total_low:5.2f} %)")

    rows = [
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r" & EUR m & share \\",
        r"\midrule",
        r"\multicolumn{3}{l}{\emph{A. The program's bill at purchasers' prices}}\\",
        rf"imported final use & ${bill['imp']:,.0f}$ & ${100*sh_imp:.1f}$ \\".replace(",", "{,}"),
        rf"product taxes & ${bill['tax']:,.0f}$ & ${100*sh_tax:.1f}$ \\".replace(",", "{,}"),
        rf"domestic production & ${bill['dom']:,.0f}$ & ${100*(1-sh_imp-sh_tax):.1f}$ \\".replace(",", "{,}"),
        rf"\emph{{program}} & ${G0:,.0f}$ & $100.0$ \\".replace(",", "{,}"),
        r"\midrule",
        r"\multicolumn{3}{l}{\emph{B. The induced gross output, by cost structure}}\\",
        rf"domestic intermediates & ${induced['dom_int']:,.0f}$ & ${100*induced['dom_int']/x.sum():.1f}$ \\".replace(",", "{,}"),
        rf"imported intermediates & ${induced['imp_int']:,.0f}$ & ${100*induced['imp_int']/x.sum():.1f}$ \\".replace(",", "{,}"),
        rf"product taxes on intermediate use & ${induced['tax_int']:,.0f}$ & ${100*induced['tax_int']/x.sum():.1f}$ \\".replace(",", "{,}"),
        rf"value added & ${induced['va']:,.0f}$ & ${100*induced['va']/x.sum():.1f}$ \\".replace(",", "{,}"),
        rf"\emph{{induced gross output}} & ${x.sum():,.0f}$ & $100.0$ \\".replace(",", "{,}"),
        r"\midrule",
        r"\multicolumn{3}{l}{\emph{C. Totals per euro of program spending}}\\",
        rf"imports, direct and embodied & ${totals['imp']:,.0f}$ & ${totals['imp']/G0:.2f}$ \\".replace(",", "{,}"),
        rf"product taxes, on the bill and on intermediates & ${totals['tax']:,.0f}$ & ${totals['tax']/G0:.2f}$ \\".replace(",", "{,}"),
        rf"value added generated & ${induced['va']:,.0f}$ & ${induced['va']/G0:.2f}$ \\".replace(",", "{,}"),
        rf"induced gross output & ${x.sum():,.0f}$ & ${x.sum()/G0:.2f}$ \\".replace(",", "{,}"),
        r"\bottomrule",
        r"\end{tabular}",
        r"\caption{The program's spending in two decompositions and a totals panel. Panel A decomposes the",
        r"program itself and sums to its size. Panel B decomposes the gross output it",
        r"calls forth and sums to that output; the split is the sectoral column identity",
        r"of the table, so it closes exactly. Panel C reports totals per euro of program",
        rf"spending, which are \emph{{not}} shares of a single base: the induced output is",
        rf"${x.sum()/G0:.2f}$ euro per euro of program, of which ${induced['va']/x.sum():.2f}$ is value added per euro of induced output and ${induced['va']/G0:.2f}$ per euro of program. The",
        r"import share of the bill in panel A is an assumption, not a datum: the published",
        r"table reports imported final use by final-demand category only, so the program",
        r"inherits the import and tax shares of the investment columns.}",
        r"\label{tab:pipeline-composition}",
        r"\end{table}",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "% GENERATED by cbase2/scripts/composition.py -- do not edit by hand.\n"
        "% Source: data/I-O_DE2019_formatiert.csv (label lookup) + data/impulses.csv (2024 row)\n"
        f"% Canary on regeneration: max |(I-P)y - (y-Zrow)| = {canary:.2e}\n"
    )
    OUT.write_text(header + "\n".join(rows) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
