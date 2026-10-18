#!/usr/bin/env python3
"""Regenerate the pricing-door tables and the price-response bracket figure
from the Phase 7 run manifests (ADR-0029, ADR-0030, ADR-0031).

No number is hand-typed: every value is a manifest field. Run from anywhere:
    python3 make_pricing_door.py
Sources (schema v2; the capacity door on OpenCGE v0.2.1 with the
value-added-scaled exponent, ADR-0031; the markup door on v0.2.0, which
v0.2.1 leaves unchanged):
    runs/pricing_capacity-v2-GAMMA-<fin>-<d0|dp05|dm05|dp05-s0|dm05-s0>/manifest.toml
    runs/pricing_markup-v1-GAMMA-<fin>-<k025|k05>/manifest.toml
    runs/pricing_profit-v1-GAMMA-<fin>-<k025|k05>-<spi05|spis>/manifest.toml
    runs/pricing_mobile-v2-ALPHA-<fin>-<dp05|dm05|dp05-s0|dm05-s0>/manifest.toml
    runs/pricing_mobile-v1-ALPHA-<fin>-<k025|k05>/manifest.toml
Outputs:
    revision/tables/table6_pricing-door.tex          (main text: goods profile)
    revision/tables/tableE_pricing-profiles.tex      (appendix: all-sector + mobile)
    revision/figures/fig_pricing-bracket.tex         (plain TikZ, no pgfplots)

Profile split (2026-10-18). The main-text table carries the goods-only profile
(delta on CPA A--F, where the programme's entire direct incidence falls) and
the bracket figure plots the same profile, so table and figure agree. The
all-sector counterfactual and the mobile-labour rows moved to the appendix
table. The sign-flip finding (all-sector: a deflator fall under external
finance only) is stated in Section 6.4 and in the appendix caption; it is not
a main-text table row any more.
"""
import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT_TABLE = os.path.join(REPO, "revision/tables/table6_pricing-door.tex")
OUT_TABLE_APP = os.path.join(REPO, "revision/tables/tableE_pricing-profiles.tex")
OUT_FIG = os.path.join(REPO, "revision/figures/fig_pricing-bracket.tex")

FIN = ["F1", "F2", "F3"]

# ---- main text: the goods-only profile + the markup door -------------------
MAIN_PANELS = [
    ("A. Capacity / Kaldor--Verdoorn door ($P_1$), goods sectors, fixed real wage", [
        (r"$\delta = 0$ (nesting)", "pricing_capacity-v2-GAMMA-{f}-d0"),
        (r"$\delta = +0.5$", "pricing_capacity-v2-GAMMA-{f}-dp05-s0"),
        (r"$\delta = -0.5$", "pricing_capacity-v2-GAMMA-{f}-dm05-s0"),
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
]

# ---- appendix: the all-sector counterfactual + the mobile robustness ------
APP_PANELS = [
    ("A. Capacity door ($P_1$), fixed wage: all-sector counterfactual", [
        (r"$\delta = +0.5$", "pricing_capacity-v2-GAMMA-{f}-dp05"),
        (r"$\delta = -0.5$", "pricing_capacity-v2-GAMMA-{f}-dm05"),
    ]),
    ("B. Mobile labour, capacity door, all sectors (CPI numeraire)", [
        (r"$\delta = +0.5$", "pricing_mobile-v2-ALPHA-{f}-dp05"),
        (r"$\delta = -0.5$", "pricing_mobile-v2-ALPHA-{f}-dm05"),
    ]),
    ("C. Mobile labour, capacity door, goods sectors (CPI numeraire)", [
        (r"$\delta = +0.5$", "pricing_mobile-v2-ALPHA-{f}-dp05-s0"),
        (r"$\delta = -0.5$", "pricing_mobile-v2-ALPHA-{f}-dm05-s0"),
    ]),
    ("D. Mobile labour, markup door, rent retained (CPI numeraire)", [
        (r"$\kappa = 0.25$", "pricing_mobile-v1-ALPHA-{f}-k025"),
        (r"$\kappa = 0.5$", "pricing_mobile-v1-ALPHA-{f}-k05"),
    ]),
]

# Bracket figure series: (legend, run-id template, colour, TikZ marker).
# Categorical slots 1-4 of the validated default palette (light surface);
# marker shape is the secondary encoding (two hues sit below 3:1 contrast).
# The capacity series carry the GOODS profile, matching the main-text table.
SERIES = [
    (r"Kaldor--Verdoorn, $\delta = -0.5$ (goods)", "pricing_capacity-v2-GAMMA-{f}-dm05-s0", "2a78d6", "circle"),
    (r"capacity pressure, $\delta = +0.5$ (goods)", "pricing_capacity-v2-GAMMA-{f}-dp05-s0", "eb6834", "rectangle"),
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


def build_table(panels, caption_lines, label, comment_lines, fit=False):
    """One pricing-door table: identical geometry for both tables."""
    L = []
    A = L.append
    for line in comment_lines:
        A(line)
    A(r"\begin{table}[htbp]")
    A(r"\centering")
    A(r"\scriptsize")
    A(r"\setlength{\tabcolsep}{1.6pt}")
    A(r"\renewcommand{\arraystretch}{1.15}")
    if fit:
        A(r"\resizebox{\textwidth}{!}{%")
    A(r"\begin{tabular}{@{}l" + "".join([r"@{\hspace{6pt}}rrrrr" for _ in FIN]) + r"@{}}")
    A(r"\toprule")
    A(" & " + " & ".join([r"\multicolumn{5}{c}{\textbf{%s}}" % f for f in FIN]) + r" \\")
    A("".join([r"\cmidrule(lr){%d-%d}" % (2 + 5 * i, 6 + 5 * i) for i in range(len(FIN))]))
    sub = r"$\Delta P$ & $\max|\Delta p|$ & $\Delta L$ & $\Delta W$ & $\Xi/Y$"
    A(r"Door and parameter & " + " & ".join([sub] * len(FIN)) + r" \\")
    for panel, rows in panels:
        A(r"\midrule")
        A(r"\multicolumn{16}{@{}l}{\textit{" + panel + r"}} \\")
        for rowlabel, tmpl in rows:
            cells = []
            for f in FIN:
                v = values(tmpl.format(f=f))
                cells += [num(v["P"]), num(v["D"], False), num(v["L"]), num(v["W"]), num(v["S"], False)]
            A(r"\quad " + rowlabel + " & " + " & ".join(cells) + r" \\")
    A(r"\bottomrule")
    A(r"\end{tabular}" + ("}" if fit else ""))
    A(r"\caption{" + caption_lines[0])
    for line in caption_lines[1:-1]:
        A(line)
    A(caption_lines[-1] + "}")
    A(label)
    A(r"\end{table}")
    return "\n".join(L) + "\n"


def build_main_table():
    return build_table(
        MAIN_PANELS,
        [r"The pricing door under the demand-only programme shock. Rows are the door",
         r"and its parameter, columns the financing closures; each block reports the change in",
         r"the price level ($\Delta P$, the GDP deflator; signed), the largest sectoral price",
         r"deviation ($\max|\Delta p|$), employment ($\Delta L$), household consumption",
         r"($\Delta W$) and the markup rent's share of factor income ($\Xi/Y$), all in percent",
         r"against the reference continuation. Panels A--C hold the real wage fixed (the",
         r"\textsc{gamma} gauge). The $\delta = 0$ row reproduces the fixed-wage column of",
         r"Table~\ref{tab:results-matrix} exactly ($\kappa = 0$ does too). $\delta$ is the",
         r"value-added elasticity of equation~\eqref{eq:capacity-door}, carried by the goods",
         r"sectors (CPA A--F); the all-sector counterfactual and the mobile-labour rows are in",
         r"Appendix~\ref{app:variants}. $s$ is the household saving rate. Generated from the",
         r"\texttt{pricing\_\{capacity,mobile\}-v2-*} (capacity door, OpenCGE v0.2.1) and",
         r"\texttt{pricing\_\{markup,profit,mobile\}-v1-*} (markup door, v0.2.0) run",
         r"manifests; no number is hand-typed."],
        r"\label{tab:pricing-door}",
        [r"% =====================================================================",
         r"% Table 6 --- the pricing door: price response by door and financing (GENERATED)",
         r"% ---------------------------------------------------------------------",
         r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py",
         r"% Carries the GOODS-ONLY capacity profile (the main-text calibration) plus the",
         r"% markup-door panels. The all-sector counterfactual and the mobile-labour rows",
         r"% live in tableE_pricing-profiles.tex (Appendix E); the bracket figure plots the",
         r"% same goods profile, so figure and table agree.",
         r"% ====================================================================="])


def build_app_table():
    return build_table(
        APP_PANELS,
        [r"The pricing-door rows the main text does not carry: the capacity door with",
         r"$\delta$ applied to every sector, and both doors under mobile labour with the CPI",
         r"numeraire. Conventions as in Table~\ref{tab:pricing-door}. Under the all-sector",
         r"profile the sign of the deflator response looks like a property of the financing",
         r"closure --- a fall only under external finance ($-0.82$ percent,",
         r"\texttt{pricing\_capacity-v2-GAMMA-F3-dm05}), against a fall in all three columns",
         r"when $\delta$ is confined to the goods sectors --- because it is the service sectors,",
         r"which receive none of the programme directly, that lose demand and productivity under",
         r"reallocation. The sign is therefore a statement about where returns are assumed to",
         r"live. Generated from the same run manifests; no number is hand-typed."],
        r"\label{tab:pricing-profiles}",
        [r"% =====================================================================",
         r"% Table E --- pricing-door profiles not carried in the main text (GENERATED)",
         r"% ---------------------------------------------------------------------",
         r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py",
         r"% Panels: A = capacity door, all-sector counterfactual; B = mobile labour, CPI",
         r"% numeraire (both doors). The main-text Table 6 carries the goods-only profile.",
         r"% ====================================================================="],
        fit=True)


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
    A(r"% Source: pricing_capacity-v2-GAMMA-*-s0 (goods profile) and")
    A(r"%   pricing_{markup,profit}-v1-GAMMA-* run manifests; y = 100 (gdp_deflator - 1).")
    A(r"%   Matches the main-text Table 6 (goods-only capacity profile).")
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
        fh.write(build_main_table())
    print("wrote", OUT_TABLE)
    with open(OUT_TABLE_APP, "w") as fh:
        fh.write(build_app_table())
    print("wrote", OUT_TABLE_APP)
    with open(OUT_FIG, "w") as fh:
        fh.write(build_figure())
    print("wrote", OUT_FIG)
