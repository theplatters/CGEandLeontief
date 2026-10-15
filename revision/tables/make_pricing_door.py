#!/usr/bin/env python3
"""Regenerate the pricing-door table and the price-response bracket figure
from the Phase 7 run manifests (ADR-0029, ADR-0030).

No number is hand-typed: every value is a manifest field. Run from anywhere:
    python3 make_pricing_door.py
Sources (schema v2, OpenCGE v0.2.0 generation):
    runs/pricing_capacity-v1-GAMMA-<fin>-<d0|dp05|dm05>/manifest.toml
    runs/pricing_markup-v1-GAMMA-<fin>-<k025|k05>/manifest.toml
    runs/pricing_profit-v1-GAMMA-<fin>-<k025|k05>-<spi05|spis>/manifest.toml
    runs/pricing_mobile-v1-ALPHA-<fin>-<dp05|dm05|k025|k05>/manifest.toml
Outputs:
    revision/tables/table6_pricing-door.tex
    revision/figures/fig_pricing-bracket.tex   (plain TikZ, no pgfplots)
"""
import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT_TABLE = os.path.join(REPO, "revision/tables/table6_pricing-door.tex")
OUT_FIG = os.path.join(REPO, "revision/figures/fig_pricing-bracket.tex")

FIN = ["F1", "F2", "F3"]

# (panel, row label, run-id template with {f} for the financing column)
PANELS = [
    ("A. Capacity / Kaldor--Verdoorn door ($P_1$), fixed real wage", [
        (r"$\delta = 0$ (nesting)", "pricing_capacity-v1-GAMMA-{f}-d0"),
        (r"$\delta = +0.5$", "pricing_capacity-v1-GAMMA-{f}-dp05"),
        (r"$\delta = -0.5$", "pricing_capacity-v1-GAMMA-{f}-dm05"),
    ]),
    ("B. Markup door ($P_2$), rent retained, fixed real wage", [
        (r"$\kappa = 0.25$", "pricing_markup-v1-GAMMA-{f}-k025"),
        (r"$\kappa = 0.5$", "pricing_markup-v1-GAMMA-{f}-k05"),
    ]),
    ("C. Markup door ($P_2$), rent distributed, fixed real wage", [
        (r"$\kappa = 0.25$, $s_\pi = 0.5$", "pricing_profit-v1-GAMMA-{f}-k025-spi05"),
        (r"$\kappa = 0.25$, $s_\pi = s$", "pricing_profit-v1-GAMMA-{f}-k025-spis"),
        (r"$\kappa = 0.5$, $s_\pi = 0.5$", "pricing_profit-v1-GAMMA-{f}-k05-spi05"),
        (r"$\kappa = 0.5$, $s_\pi = s$", "pricing_profit-v1-GAMMA-{f}-k05-spis"),
    ]),
    ("D. Robustness: mobile labour, CPI numeraire", [
        (r"$\delta = +0.5$", "pricing_mobile-v1-ALPHA-{f}-dp05"),
        (r"$\delta = -0.5$", "pricing_mobile-v1-ALPHA-{f}-dm05"),
        (r"$\kappa = 0.25$, retained", "pricing_mobile-v1-ALPHA-{f}-k025"),
        (r"$\kappa = 0.5$, retained", "pricing_mobile-v1-ALPHA-{f}-k05"),
    ]),
]

# Bracket figure series: (legend, run-id template, colour, TikZ marker).
# Categorical slots 1-4 of the validated default palette (light surface);
# marker shape is the secondary encoding (two hues sit below 3:1 contrast).
SERIES = [
    (r"Kaldor--Verdoorn, $\delta = -0.5$", "pricing_capacity-v1-GAMMA-{f}-dm05", "2a78d6", "circle"),
    (r"capacity pressure, $\delta = +0.5$", "pricing_capacity-v1-GAMMA-{f}-dp05", "eb6834", "rectangle"),
    (r"markup retained, $\kappa = 0.5$", "pricing_markup-v1-GAMMA-{f}-k05", "1baf7a", "triangle"),
    (r"markup distributed, $\kappa = 0.5$, $s_\pi = s$", "pricing_profit-v1-GAMMA-{f}-k05-spis", "eda100", "diamond"),
]


def manifest(rid):
    with open(f"{REPO}/runs/{rid}/manifest.toml", "rb") as fh:
        d = tomllib.load(fh)
    assert d["status"] == "executed" and d["gates"]["overall"] == "pass", rid
    return d


def values(rid):
    m = manifest(rid)["metrics"]
    return dict(P=100 * (m["gdp_deflator"] - 1.0), D=100 * m["max_abs_price_dev"],
                L=100 * (m["employment"] - 1.0), W=100 * m["consumption_rel"],
                S=100 * m.get("profit_share", 0.0))


def num(x, signed=True):
    if abs(x) < 5e-3:
        return r"$0.00$"
    if not signed:
        return "$%.2f$" % x
    return ("$+%.2f$" % x) if x > 0 else ("$%.2f$" % x)


def build_table():
    L = []
    A = L.append
    A(r"% =====================================================================")
    A(r"% Table 6 --- the pricing door: price response by door and financing (GENERATED)")
    A(r"% ---------------------------------------------------------------------")
    A(r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py")
    A(r"% Source: runs/pricing_{capacity,markup,profit,mobile}-v1-*/manifest.toml (schema v2,")
    A(r"%   OpenCGE v0.2.0). Per financing column: \Delta P GDP deflator, max|p-1| price")
    A(r"%   dispersion, \Delta L employment, \Delta W household consumption, \Pi/Y profit share;")
    A(r"%   all in percent vs. the reference continuation.")
    A(r"% =====================================================================")
    A(r"\begin{table}[htbp]")
    A(r"\centering")
    A(r"\scriptsize")
    A(r"\setlength{\tabcolsep}{1.6pt}")
    A(r"\renewcommand{\arraystretch}{1.15}")
    A(r"\begin{tabular}{@{}l" + "".join([r"@{\hspace{6pt}}rrrrr" for _ in FIN]) + r"@{}}")
    A(r"\toprule")
    A(" & " + " & ".join([r"\multicolumn{5}{c}{\textbf{%s}}" % f for f in FIN]) + r" \\")
    A("".join([r"\cmidrule(lr){%d-%d}" % (2 + 5 * i, 6 + 5 * i) for i in range(len(FIN))]))
    sub = r"$\Delta P$ & $\max|\Delta p|$ & $\Delta L$ & $\Delta W$ & $\Pi/Y$"
    A(r"Door and parameter & " + " & ".join([sub] * len(FIN)) + r" \\")
    for panel, rows in PANELS:
        A(r"\midrule")
        A(r"\multicolumn{16}{@{}l}{\textit{" + panel + r"}} \\")
        for label, tmpl in rows:
            cells = []
            for f in FIN:
                v = values(tmpl.format(f=f))
                cells += [num(v["P"]), num(v["D"], False), num(v["L"]), num(v["W"]), num(v["S"], False)]
            A(r"\quad " + label + " & " + " & ".join(cells) + r" \\")
    A(r"\bottomrule")
    A(r"\end{tabular}")
    A(r"\caption{The pricing door under the demand-only programme shock. Rows are the door")
    A(r"and its parameter, columns the financing closures; each block reports the change in")
    A(r"the price level ($\Delta P$, the GDP deflator; signed), the largest sectoral price")
    A(r"deviation ($\max|\Delta p|$), employment ($\Delta L$), household consumption")
    A(r"($\Delta W$) and the markup rent's share of factor income ($\Pi/Y$), all in percent")
    A(r"against the reference continuation. Panels A--C hold the real wage fixed (the")
    A(r"\textsc{gamma} gauge); panel D is the mobile-labour robustness row under the CPI")
    A(r"numeraire, where employment is fixed and prices move relative to the CPI. The")
    A(r"$\delta = 0$ row reproduces the fixed-wage column of Table~\ref{tab:results-matrix}")
    A(r"exactly ($\kappa = 0$ does too). $s$ is the household saving rate. Generated from the")
    A(r"\texttt{pricing\_\{capacity,markup,profit,mobile\}-v1-*} run manifests (OpenCGE")
    A(r"v0.2.0); no number is hand-typed.}")
    A(r"\label{tab:pricing-door}")
    A(r"\end{table}")
    return "\n".join(L) + "\n"


def marker(mk, size):
    """TikZ node style for a marker; regular polygons size by their
    circumcircle, so the triangle and diamond get a larger box to match the
    circle's visual weight."""
    shape = {"triangle": "regular polygon, regular polygon sides=3",
             "diamond": "diamond"}.get(mk, mk)
    scale = {"triangle": 1.45, "diamond": 1.3}.get(mk, 1.0)
    return "%s, minimum size=%.2fmm" % (shape, size * scale)


def build_figure():
    # Dot plot: x = financing column, y = GDP-deflator change in percent.
    data = [[values(t.format(f=f))["P"] for f in FIN] for _, t, _, _ in SERIES]
    lo = min(min(r) for r in data)
    hi = max(max(r) for r in data)
    ymin, ymax = (int(lo // 0.5) * 0.5) - 0.5, (int(hi // 0.5) + 1) * 0.5
    H = 6.0                      # plot height in cm
    sy = H / (ymax - ymin)
    y = lambda v: (v - ymin) * sy
    xs = [1.5, 4.5, 7.5]         # column centres in cm
    dx = [-0.6, -0.2, 0.2, 0.6]  # per-series offsets
    L = []
    A = L.append
    A(r"% =====================================================================")
    A(r"% Figure --- the price-response bracket (GENERATED)")
    A(r"% Do NOT edit by hand. Regenerate with:  revision/tables/make_pricing_door.py")
    A(r"% Source: the pricing_*-v1-GAMMA-* run manifests (OpenCGE v0.2.0); y = 100 (gdp_deflator - 1).")
    A(r"% Colours: validated categorical slots 1-4 (light); markers are the secondary encoding.")
    A(r"% =====================================================================")
    for i, (_, _, col, _) in enumerate(SERIES):
        A(r"\definecolor{pd%d}{HTML}{%s}" % (i + 1, col))
    A(r"\begin{tikzpicture}[font=\footnotesize]")
    # recessive grid and axis
    t = ymin
    while t <= ymax + 1e-9:
        A(r"\draw[black!12, line width=0.4pt] (0,%.3f) -- (9,%.3f);" % (y(t), y(t)))
        A(r"\node[anchor=east, text=black!60] at (0,%.3f) {$%+.1f$};" % (y(t), t) if abs(t) > 1e-9
          else r"\node[anchor=east, text=black!60] at (0,%.3f) {$0.0$};" % y(t))
        t += 0.5
    A(r"\draw[black!55, line width=0.6pt] (0,%.3f) -- (9,%.3f);" % (y(0), y(0)))
    A(r"\node[rotate=90, anchor=south, text=black!70] at (-0.9,%.3f) {GDP deflator, \%% change};" % (H / 2))
    for j, f in enumerate(FIN):
        A(r"\node[anchor=north, text=black!70] at (%.2f,-0.1) {%s};" % (xs[j], f))
        col_vals = [data[i][j] for i in range(len(SERIES))]
        # the bracket: the span of the four doors in this financing column
        A(r"\draw[black!35, line width=1.2pt] (%.2f,%.3f) -- (%.2f,%.3f);"
          % (xs[j] - 0.95, y(min(col_vals)), xs[j] - 0.95, y(max(col_vals))))
        for i, (_, _, _, mk) in enumerate(SERIES):
            A(r"\node[draw=white, line width=0.6pt, fill=pd%d, %s, inner sep=0pt] at (%.2f,%.3f) {};"
              % (i + 1, marker(mk, 2.4), xs[j] + dx[i], y(data[i][j])))
    # legend (always present: four series, identity never colour alone)
    for i, (lab, _, _, mk) in enumerate(SERIES):
        yy = H - 0.05 - 0.42 * i
        A(r"\node[draw=white, line width=0.6pt, fill=pd%d, %s, inner sep=0pt] at (9.5,%.3f) {};"
          % (i + 1, marker(mk, 2.2), yy))
        A(r"\node[anchor=west, text=black!80] at (9.7,%.3f) {%s};" % (yy, lab))
    A(r"\node[anchor=west, text=black!60] at (9.4,%.3f) {grey bar: the bracket across doors};" % (H - 0.05 - 0.42 * 4))
    A(r"\end{tikzpicture}")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    with open(OUT_TABLE, "w") as fh:
        fh.write(build_table())
    print("wrote", OUT_TABLE)
    with open(OUT_FIG, "w") as fh:
        fh.write(build_figure())
    print("wrote", OUT_FIG)
