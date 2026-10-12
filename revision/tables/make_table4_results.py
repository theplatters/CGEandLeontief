#!/usr/bin/env python3
"""Regenerate table4_results-5x3.tex from the v10 run manifests.

No number is hand-typed: every value is a manifest field. Run from anywhere:
    python3 make_table4_results.py
Source: runs/matrix_5x3-v10-<closure>-<fin>[-<variant>]/manifest.toml (schema v2).
Columns are in table3c order (no-reallocation leftmost). The elastic-supply
column reads the SECTORAL cells (eta_{s,i} = 0.5, variant ``etas05``; ADR-0022),
not the single-wage BETA-F* benchmark, which coincides with the standard CGE
column under demand-only shocks (Appendix B, Prop. ab).
"""
import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT  = os.path.join(REPO, "revision/tables/table4_results-5x3.tex")

# run-id variant suffix per column ("" = the plain matrix cell)
variant = {"BETA": "-etas05"}

# column order matches table3c (no-reallocation leftmost)
cols = [("BF",    "No-reallocation",      "immobile labour, sectoral wages"),
        ("ALPHA", "Standard CGE",         "one flexible wage, mobile labour"),
        ("BETA",  "Elastic labour supply","sectoral supply on the real wage ($\\eta_{s,i} = 0.5$)"),
        ("GAMMA", "Fixed real wage",      "real wage pegged, employment uncapped"),
        ("DELTA", "Leontief closure",     "real wage pegged, quantities absorb")]
fin = [("F1", "repurposing of consumption"),
       ("F2", "refinancing via taxes"),
       ("F3", "external / debt finance")]


def get(c, f):
    d = tomllib.load(open(f"{REPO}/runs/matrix_5x3-v10-{c}-{f}{variant.get(c, '')}/manifest.toml", "rb"))
    m = d["metrics"]
    return dict(Y=100 * m["gdp_rel"], L=100 * (m["employment"] - 1.0),
                P=100 * (m["gdp_deflator"] - 1.0),
                D=100 * m["external_position"], W=100 * m["consumption_rel"])


def num(x):
    if abs(x) < 5e-5:
        return r"$0.00$"
    return ("$+%.2f$" % x) if x > 0 else ("$%.2f$" % x)


def cell(v):
    return (r"\begin{tabular}{@{}l@{\hspace{3pt}}r@{}}" "\n"
            r"$\Delta Y$ & " + num(v["Y"]) + r" \\" "\n"
            r"$\Delta L$ & " + num(v["L"]) + r" \\" "\n"
            r"$\Delta P$ & " + num(v["P"]) + r" \\" "\n"
            r"$\Delta D$ & " + num(v["D"]) + r" \\" "\n"
            r"$\Delta W$ & " + num(v["W"]) + r"\end{tabular}")


def build():
    L = []
    A = L.append
    A(r"% =====================================================================")
    A(r"% Table 4 --- results across the 5x3 closure matrix (GENERATED)")
    A(r"% ---------------------------------------------------------------------")
    A(r"% Do NOT edit by hand. Regenerate with:  make_table4_results.py")
    A(r"% Source: runs/matrix_5x3-v10-<closure>-<fin>/manifest.toml (schema v2);")
    A(r"%   elastic-supply column: the sectoral cells matrix_5x3-v10-BETA-<fin>-etas05.")
    A(r"% Columns in table3c order (no-reallocation leftmost). Five outcomes per cell:")
    A(r"%   \Delta Y real GDP, \Delta L employment, \Delta P price level (GDP deflator),")
    A(r"%   \Delta D net external position (= F + B_gov, 'external debt'),")
    A(r"%   \Delta W household consumption (welfare).")
    A(r"% All in percent vs. the reference continuation; GDP at basic prices = 1.")
    A(r"% =====================================================================")
    A(r"\begin{table}[htbp]")
    A(r"\centering")
    A(r"\footnotesize")
    A(r"\setlength{\tabcolsep}{3pt}")
    A(r"\setlength{\parindent}{0pt}")
    A(r"\setlength{\parskip}{2pt}")
    A(r"\renewcommand{\arraystretch}{1.25}")
    A(r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{1.7cm}"
      + "".join([r" >{\raggedright\arraybackslash}p{2.5cm}" for _ in cols]) + r"@{}}")
    A(r"\toprule")
    A(r"\textbf{Financing closure}")
    for _, name, desc in cols:
        A(r"& \textbf{" + name + r"}\newline \scriptsize " + desc)
    A(r"\\")
    A(r"\midrule")
    for i, (f, fd) in enumerate(fin):
        A(r"\textbf{" + f + r"}\newline \scriptsize " + fd)
        for c, _, _ in cols:
            A(r"& " + cell(get(c, f)))
        A(r"\\")
        if i < len(fin) - 1:
            A(r"\midrule")
    A(r"\bottomrule")
    A(r"\end{tabular}")
    A(r"\caption{Results across the closure matrix: the three financing closures (rows)")
    A(r"crossed with the five labour-market closures (columns). Each cell reports the")
    A(r"change in real GDP ($\Delta Y$), employment ($\Delta L$), the price level")
    A(r"($\Delta P$, the GDP deflator), the net external position")
    A(r"($\Delta D = F + B_{\mathrm{gov}}$, the external-account block; positive = net external")
    A(r"financing received, i.e.\ external debt) and household consumption welfare")
    A(r"($\Delta W$), all in percent against the reference continuation. Generated from the")
    A(r"\texttt{matrix\_5x3-v10-*} run manifests; no number is hand-typed. The")
    A(r"elastic-supply column is the sectoral closure at $\eta_{s,i} = 0.5$")
    A(r"(\texttt{matrix\_5x3-v10-BETA-F*-etas05}). The fixed-wage and Leontief columns")
    A(r"coincide (bit-identical), and F2 and F3 coincide in the three flexible-wage")
    A(r"columns, so the fifteen cells describe nine distinct economies: the standard-CGE")
    A(r"column is exactly output-neutral, and equilibrium prices move away from the")
    A(r"baseline only in the two sectoral-wage columns (no reallocation, elastic supply).}")
    A(r"\label{tab:results-matrix}")
    A(r"\end{table}")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    open(OUT, "w").write(build())
    print("wrote", OUT)
