#!/usr/bin/env python3
"""Regenerate table4_results-5x3.tex from the v10 run manifests.

No number is hand-typed: every value is a manifest field. Run from anywhere:
    python3 make_table4_results.py
Source: runs/matrix_5x3-v10-<closure>-<fin>[-<variant>]/manifest.toml (schema v2).
Columns are in table3c order (no-reallocation leftmost). The elastic-supply
column reads the SECTORAL cells at the four executed uniform elasticities
eta_{s,i} in {0.25, 0.5, 1, 2} (variants ``etas025`` ... ``etas2``; ADR-0022;
0.5 is the headline value), one sub-column each, not the single-wage BETA-F*
benchmark, which coincides with the standard CGE
column under demand-only shocks (Appendix C, Prop. ab).
"""
import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT  = os.path.join(REPO, "revision/tables/table4_results-5x3.tex")

# run-id variant suffix per column ("" = the plain matrix cell)
variant = {"BETA": "-etas05"}

# elastic-supply sub-columns: (label, run-id suffix); the headline value is bold
etas = [("0.25", "-etas025"), ("0.5", "-etas05"), ("1", "-etas1"), ("2", "-etas2")]
ETAS_HEAD = "0.5"

# column order matches table3c (no-reallocation leftmost)
cols = [("BF",    "No-reallocation",      "immobile labour, sectoral wages"),
        ("ALPHA", "Standard CGE",         "one flexible wage, mobile labour"),
        ("BETA",  "Elastic labour supply","sectoral supply on the real wage, by $\\eta_{s,i}$"),
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


OUTCOMES = [("Y", r"$\Delta Y$"), ("L", r"$\Delta L$"), ("P", r"$\Delta P$"),
            ("D", r"$\Delta D$"), ("W", r"$\Delta W$")]
# inner layout of the elastic-supply column: outcome label + one sub-column per eta
# (0.6 + 4 x (0.105 + 0.76) = 4.06cm, inside the 4.2cm column: a gap before the next one)
ETAS_SPEC = (r"@{}>{\raggedright\arraybackslash}p{0.6cm}"
             r"*{%d}{@{\hspace{3pt}}>{\raggedleft\arraybackslash}p{0.76cm}}@{}" % len(etas))


def etas_head():
    lab = [(r"$\mathbf{%s}$" if e == ETAS_HEAD else r"$%s$") % e for e, _ in etas]
    # footnotesize, as the cells below, so the sub-columns line up
    return (r"{\footnotesize\begin{tabular}{" + ETAS_SPEC + "}\n"
            r"$\eta_{s,i}$ & " + " & ".join(lab) + r"\end{tabular}}")


def etas_cell(f):
    vs = [get("BETA", f, suf) for _, suf in etas]
    rows = [lab + " & " + " & ".join(num(v[k]) for v in vs) for k, lab in OUTCOMES]
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
    A(r"%   elastic-supply column: the sectoral cells matrix_5x3-v10-BETA-<fin>-etas{025,05,1,2},")
    A(r"%   one sub-column per uniform eta_{s,i} (0.5 = headline value, bold).")
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
    width = {"BETA": "4.2cm"}
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
    A(r"($\Delta D = R + B_{\mathrm{gov}}$, the external-account block; positive = net external")
    A(r"financing received, i.e.\ external debt) and household consumption welfare")
    A(r"($\Delta W$), all in percent against the reference continuation. Generated from the")
    A(r"\texttt{matrix\_5x3-v10-*} run manifests; no number is hand-typed. The")
    A(r"elastic-supply column is the sectoral closure at four uniform supply elasticities,")
    A(r"$\eta_{s,i} \in \{0.25, 0.5, 1, 2\}$ (\texttt{matrix\_5x3-v10-BETA-F*-etas025},")
    A(r"\dots, \texttt{-etas2}); $\eta_{s,i} = 0.5$ (bold) is the headline value used in")
    A(r"the text. The fixed-wage and Leontief columns")
    A(r"coincide (bit-identical), and F2 and F3 coincide in the three flexible-wage")
    A(r"columns, so at the headline elasticity the fifteen cells describe nine distinct economies: the standard-CGE")
    A(r"column is exactly output-neutral, and equilibrium prices move away from the")
    A(r"baseline only in the two sectoral-wage columns (no reallocation, elastic supply).}")
    A(r"\label{tab:results-matrix}")
    A(r"\end{table}")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    open(OUT, "w").write(build())
    print("wrote", OUT)
