#!/usr/bin/env python3
"""Regenerate table4_results-5x3.tex from the v10 run manifests.

No number is hand-typed: every value is a manifest field. Run from anywhere:
    python3 make_table4_results.py
Source: runs/matrix_5x3-v10-<closure>-<fin>[-<variant>]/manifest.toml (schema v2).
Columns are in table3c order (no-reallocation leftmost).

The elastic-supply column reports a BRACKET rather than one sub-column per
executed elasticity: each outcome is the min--max interval over
eta_{s,i} in ETA_BRACKET, read from the executed sectoral cells
``-etas05`` ... ``-etas2`` (ADR-0022). The headline point value is the lower
end of the bracket, eta_{s,i} = 0.5.
"""

import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT  = os.path.join(REPO, "revision/tables/table4_results-5x3.tex")

# run-id variant suffix per column ("" = the plain matrix cell)
variant = {"BETA": "-etas05"}

# the bracket: executed sectoral cells, low end first (0.5 = headline value)
ETA_BRACKET = [("0.5", "-etas05"), ("1", "-etas1"), ("2", "-etas2")]
ETA_LO, ETA_HI = ETA_BRACKET[0][0], ETA_BRACKET[-1][0]

# column order matches table3c (no-reallocation leftmost)
cols = [("BF",    "No-reallocation",      "immobile labour, sectoral wages"),
        ("ALPHA", "Standard CGE",         "one flexible wage, mobile labour"),
        ("BETA",  "Elastic labour supply","sectoral supply on the real wage"),
        ("GAMMA", "Fixed real wage",      "real wage pegged, employment uncapped"),
        ("DELTA", "Leontief closure",     "real wage pegged, quantities absorb")]
fin = [("F1", "repurposing of consumption"),
       ("F2", "refinancing via taxes"),
       ("F3", "external / debt finance")]


def get(c, f, suffix=None):
    suffix = variant.get(c, "") if suffix is None else suffix
    d = tomllib.load(open(f"{REPO}/runs/matrix_5x3-v10-{c}-{f}{suffix}/manifest.toml", "rb"))
    m = d["metrics"]
    return dict(Y=100 * m["gdp_rel"], L=100 * (m["employment"] - 1.0),
                P=100 * (m["gdp_deflator"] - 1.0),
                D=100 * m["external_position"], WC=100 * m["consumption_rel"])


def num(x):
    if abs(x) < 5e-5:
        return r"$0.00$"
    return ("$+%.2f$" % x) if x > 0 else ("$%.2f$" % x)


def num_range(lo, hi):
    """min--max interval; collapses to a point if the outcome is flat."""
    if abs(hi - lo) < 5e-5:
        return num(lo)
    return num(lo) + r"--" + num(hi)


def cell(v):
    return (r"\begin{tabular}{@{}l@{\hspace{3pt}}r@{}}" "\n"
            r"$\Delta Y$ & " + num(v["Y"]) + r" \\" "\n"
            r"$\Delta L$ & " + num(v["L"]) + r" \\" "\n"
            r"$\Delta P$ & " + num(v["P"]) + r" \\" "\n"
            r"$\Delta D$ & " + num(v["D"]) + r" \\" "\n"
            r"$\Delta W_C$ & " + num(v["WC"]) + r"\end{tabular}")


OUTCOMES = [("Y", r"$\Delta Y$"), ("L", r"$\Delta L$"), ("P", r"$\Delta P$"),
            ("D", r"$\Delta D$"), ("WC", r"$\Delta W_C$")]
# inner layout of the elastic-supply column: outcome label + one bracketed column
# (0.6 + 0.105 + 2.0 = 2.7cm inside the 2.8cm column)
ETAS_SPEC = r"@{}>{\raggedright\arraybackslash}p{0.85cm}@{\hspace{3pt}}>{\raggedleft\arraybackslash}p{2.2cm}@{}"


def etas_head():
    return (r"{\footnotesize\begin{tabular}{" + ETAS_SPEC + "}\n"
            r"& $\eta_{s,i} \in [" + ETA_LO + r"," + ETA_HI + r"]$\end{tabular}}")


def etas_cell(f):
    vs = [get("BETA", f, suf) for _, suf in ETA_BRACKET]
    rows = [lab + " & " + num_range(min(v[k] for v in vs), max(v[k] for v in vs))
            for k, lab in OUTCOMES]
    return (r"\begin{tabular}{" + ETAS_SPEC + "}\n"
            + " \\\\\n".join(rows) + r"\end{tabular}")


def build():
    L = []
    A = L.append
    A(r"% =====================================================================")
    A(r"% Table 4 --- results across the 5x3 closure matrix (GENERATED)")
    A(r"% ---------------------------------------------------------------------")
    A(r"% Do NOT edit by hand. Regenerate with:  make_table4_results.py")
    A(r"% Source: runs/matrix_5x3-v10-<closure>-<fin>/manifest.toml (schema v2);")
    A(r"%   elastic-supply column: the executed sectoral cells")
    A(r"%   matrix_5x3-v10-BETA-<fin>-etas{05,1,2}, reported as a min--max bracket.")
    A(r"% Columns in table3c order (no-reallocation leftmost). Five outcomes per cell:")
    A(r"%   \Delta Y real GDP, \Delta L employment, \Delta P price level (GDP deflator),")
    A(r"%   \Delta D net external position (= F + B_gov, 'external debt'),")
    A(r"%   \Delta W_C household consumption welfare.")
    A(r"% All in percent vs. the reference continuation; GDP at basic prices = 1.")
    A(r"% =====================================================================")
    A(r"\begin{table}[htbp]")
    A(r"\centering")
    A(r"\footnotesize")
    A(r"\setlength{\tabcolsep}{3pt}")
    A(r"\setlength{\parindent}{0pt}")
    A(r"\setlength{\parskip}{2pt}")
    A(r"\renewcommand{\arraystretch}{1.25}")
    width = {"BETA": "3.3cm"}
    A(r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{1.5cm}"
      + "".join([r" >{\raggedright\arraybackslash}p{%s}" % width.get(c, "2.1cm") for c, _, _ in cols])
      + r"@{}}")
    A(r"\toprule")
    A(r"\textbf{Financing closure}")
    for c, name, desc in cols:
        A(r"& \textbf{" + name + r"}\newline \scriptsize " + desc
          + ((r"\newline " + etas_head()) if c == "BETA" else ""))
    A(r"\\")
    A(r"\midrule")
    for i, (f, fd) in enumerate(fin):
        A(r"\textbf{" + f + r"}\newline \scriptsize " + fd)
        for c, _, _ in cols:
            A(r"& " + (etas_cell(f) if c == "BETA" else cell(get(c, f))))
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
    A(r"($\Delta W_C$), all in percent against the reference continuation. Generated from the")
    A(r"\texttt{matrix\_5x3-v10-*} run manifests; no number is hand-typed. The")
    A(r"elastic-supply column is bracketed rather than tabulated at one elasticity:")
    A(r"each entry is the min--max interval over the executed supply elasticities")
    cells = (r"\texttt{matrix\_5x3-v10-BETA-F*" + ETA_BRACKET[0][1] + "}"
             + "".join(", " + r"\texttt{" + suf + "}" for _, suf in ETA_BRACKET[1:]))
    A(r"$\eta_{s,i} \in \{" + r",\ ".join(e for e, _ in ETA_BRACKET) + r"\}$"
      + r" (" + cells + r"),")
    A(r"with $\eta_{s,i} = " + ETA_LO + r"$ the headline point value. The fixed-wage and")
    A(r"Leontief columns coincide (bit-identical), and F2 and F3 coincide in the three")
    A(r"flexible-wage columns, so the fifteen cells describe nine distinct economies:")
    A(r"the standard-CGE column is exactly output-neutral, and equilibrium prices move")
    A(r"away from the baseline only in the two sectoral-wage columns (no reallocation,")
    A(r"elastic supply).}")
    A(r"\label{tab:results-matrix}")
    A(r"\end{table}")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    open(OUT, "w").write(build())
    print("wrote", OUT)
