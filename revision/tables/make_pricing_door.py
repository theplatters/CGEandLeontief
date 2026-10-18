#!/usr/bin/env python3
"""Regenerate the pricing-door tables from the Phase 7 run manifests
(ADR-0029, ADR-0030, ADR-0031, ADR-0032).

No number is hand-typed: every value is a manifest field. Run from anywhere:
    python3 make_pricing_door.py
Sources (schema v2, all on OpenCGE v0.2.1):
    capacity door (P1), by labour row:
      BF, BETA, DELTA  runs/pricing_capacity_rows-v1-<row>-<fin>-<d0|dp05|dm05|dp05-s0|dm05-s0>
      ALPHA            runs/pricing_capacity_rows-v1-ALPHA-<fin>-d0,
                       runs/pricing_mobile-v2-ALPHA-<fin>-<dp05|dm05|dp05-s0|dm05-s0>
      GAMMA            runs/pricing_capacity-v2-GAMMA-<fin>-<d0|dp05|dm05|dp05-s0|dm05-s0>
    markup door (P2), every labour row:
      runs/pricing_markup-v2-<row>-<fin>-<k025|k05>
      runs/pricing_profit-v2-<row>-<fin>-<k025|k05>-<spi05|spis>
Outputs:
    revision/tables/table6_pricing-capacity.tex   (main text: capacity door, goods profile)
    revision/tables/table7_pricing-markup.tex     (main text: markup door)
    revision/tables/tableE_pricing-profiles.tex   (appendix: all-sector capacity, kappa = 0.25 distributed)

Split (2026-10-18, ADR-0032). The two doors get one table each, and each table
carries all five labour rows, so the dissolution of the demand-only
equivalences (fixed wage vs Leontief, standard CGE vs elastic supply) is
read off the tables rather than stated. The capacity table carries the
goods-only profile (delta on CPA A--F, where the programme's entire direct
incidence falls); the all-sector counterfactual is in the appendix table. The
price-response bracket figure was removed from the manuscript the same day.
"""
import tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT_CAP = os.path.join(REPO, "revision/tables/table6_pricing-capacity.tex")
OUT_MK = os.path.join(REPO, "revision/tables/table7_pricing-markup.tex")
OUT_TABLE_APP = os.path.join(REPO, "revision/tables/tableE_pricing-profiles.tex")

FIN = ["F1", "F2", "F3"]

# Labour rows in the order and with the names of the results matrix.
ROWS = [
    ("BF", "No-reallocation"),
    ("ALPHA", "Standard CGE"),
    ("BETA", r"Elastic labour supply ($\eta_{s,i} = 0.5$)"),
    ("GAMMA", "Fixed real wage"),
    ("DELTA", "Leontief closure"),
]


def cap_id(row, v):
    """Run-id template of the capacity cell `v` (d0, dp05, ...) on a labour row."""
    if row == "BF" and v == "dm05":
        # BF-F1-dm05 failed on a batch stall and was rerun as -v2 (ADR-0004);
        # F2/F3 ran clean under the original ids.
        return "pricing_capacity_rows-v1-BF-{f}-dm05" + "{v2}"
    if row == "GAMMA":
        return "pricing_capacity-v2-GAMMA-{f}-" + v
    if row == "ALPHA" and v != "d0":
        return "pricing_mobile-v2-ALPHA-{f}-" + v
    return "pricing_capacity_rows-v1-" + row + "-{f}-" + v


def per_row(rows_of):
    return [(name, rows_of(row)) for row, name in ROWS]


MK_PANELS = per_row(lambda row: [
    (r"retained, $\kappa = 0.25$", "pricing_markup-v2-" + row + "-{f}-k025"),
    (r"retained, $\kappa = 0.5$", "pricing_markup-v2-" + row + "-{f}-k05"),
    (r"distributed, $\kappa = 0.5$, $s_\pi = 0.5$", "pricing_profit-v2-" + row + "-{f}-k05-spi05"),
    (r"distributed, $\kappa = 0.5$, $s_\pi = s$", "pricing_profit-v2-" + row + "-{f}-k05-spis"),
])

APP_PANELS = (
    [("A. Capacity door, all sectors: " + name, [
        (r"$\delta = +0.5$", cap_id(row, "dp05")),
        (r"$\delta = -0.5$", cap_id(row, "dm05")),
    ]) for row, name in ROWS]
    + [("B. Markup door, rent distributed, $\\kappa = 0.25$: " + name, [
        (r"$s_\pi = 0.5$", "pricing_profit-v2-" + row + "-{f}-k025-spi05"),
        (r"$s_\pi = s$", "pricing_profit-v2-" + row + "-{f}-k025-spis"),
    ]) for row, name in ROWS]
)

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


COLS_CAP = [("P", r"$\Delta P$", True), ("D", r"$\max|\Delta p|$", False),
            ("L", r"$\Delta L$", True), ("W", r"$\Delta W$", True)]
COLS_MK = COLS_CAP + [("S", r"$\Xi/Y$", False)]


def build_table(panels, cols, caption_lines, label, comment_lines, fit=False):
    """One pricing-door table: a panel per labour row, a block per financing column."""
    k = len(cols)
    ncol = 1 + k * len(FIN)
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
    A(r"\begin{tabular}{@{}l" + "".join([r"@{\hspace{6pt}}" + "r" * k for _ in FIN]) + r"@{}}")
    A(r"\toprule")
    A(" & " + " & ".join([r"\multicolumn{%d}{c}{\textbf{%s}}" % (k, f) for f in FIN]) + r" \\")
    A("".join([r"\cmidrule(lr){%d-%d}" % (2 + k * i, 1 + k * (i + 1)) for i in range(len(FIN))]))
    sub = " & ".join(h for _, h, _ in cols)
    A(r"Labour closure and door & " + " & ".join([sub] * len(FIN)) + r" \\")
    for panel, rows in panels:
        A(r"\midrule")
        A(r"\multicolumn{%d}{@{}l}{\textit{%s}} \\" % (ncol, panel))
        for rowlabel, tmpl in rows:
            cells = []
            for f in FIN:
                v = values(tmpl.format(f=f, v2="-v2" if f == "F1" else ""))
                cells += [num(v[key], signed) for key, _, signed in cols]
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


FIN_NOTE = {"F1": "repurposing of consumption", "F2": "refinancing via taxes",
            "F3": "external / debt finance"}


def build_guided_table(panels, cols, stub, caption_lines, label, comment_lines):
    """The pricing-door layout tuned for reading across many columns: financing
    blocks named as in Table 4 and set apart by wider gaps, labour panels separated
    by vertical space instead of rules, and the delta = 0 reference row of every
    panel shaded as a guide line across the row (it is also the row that reproduces
    Table 4); thin grey rules between the columns of a block keep each value
    under its header (they stop at the panel headings and block gaps). The block gaps are empty spacer columns rather than @{} glue, and
    the outer columns keep their padding, so the shading runs unbroken from rule
    end to rule end."""
    k = len(cols)
    ncol = 1 + (k + 1) * len(FIN)
    L = []
    A = L.append
    for line in comment_lines:
        A(line)
    A(r"\begin{table}[htbp]")
    A(r"\centering")
    A(r"\footnotesize")
    A(r"\definecolor{tableguide}{gray}{0.92}")
    A(r"\definecolor{tableline}{gray}{0.72}")
    A(r"\setlength{\tabcolsep}{3pt}")
    A(r"\renewcommand{\arraystretch}{1.2}")
    guide = r"!{\color{tableline}\vrule width 0.3pt}"
    A(r"\begin{tabular}{l" + "".join([r"@{\hspace{3pt}}p{2pt}@{\hspace{3pt}}" + guide.join("r" * k)
                                      for _ in FIN]) + r"}")
    A(r"\toprule")
    A(" & " + " & ".join([r"& \multicolumn{%d}{c}{\textbf{%s}}" % (k, f) for f in FIN]) + r" \\")
    A(" & " + " & ".join([r"& \multicolumn{%d}{c}{\scriptsize %s}" % (k, FIN_NOTE[f]) for f in FIN])
      + r" \\")
    A("".join([r"\cmidrule(lr){%d-%d}" % (3 + (k + 1) * i, 2 + (k + 1) * i + k)
               for i in range(len(FIN))]))
    sub = " & ".join(h for _, h, _ in cols)
    A(stub + " & " + " & ".join(["& " + sub] * len(FIN)) + r" \\")
    A(r"\midrule")
    for i, (panel, rows) in enumerate(panels):
        if i:
            A(r"\addlinespace[6pt]")
        A(r"\multicolumn{%d}{l}{\textbf{%s}} \\" % (ncol, panel))
        for j, (rowlabel, tmpl) in enumerate(rows):
            cells = []
            for f in FIN:
                v = values(tmpl.format(f=f, v2="-v2" if f == "F1" else ""))
                cells += [""] + [num(v[key], signed) for key, _, signed in cols]
            shade = r"\rowcolor{tableguide}" if j == 0 else ""
            A(shade + r"\quad " + rowlabel + " & " + " & ".join(cells) + r" \\")
    A(r"\bottomrule")
    A(r"\end{tabular}")
    A(r"\caption{" + caption_lines[0])
    for line in caption_lines[1:-1]:
        A(line)
    A(caption_lines[-1] + "}")
    A(label)
    A(r"\end{table}")
    return "\n".join(L) + "\n"


CAP_GUIDED_PANELS = per_row(lambda row: [
    (r"$0$ (no door)", cap_id(row, "d0")),
    (r"$+0.5$", cap_id(row, "dp05-s0")),
    (r"$-0.5$", cap_id(row, "dm05-s0")),
])


def build_cap_table():
    return build_guided_table(
        CAP_GUIDED_PANELS, COLS_CAP, r"$\delta$",
        [r"The capacity / Kaldor--Verdoorn door under the demand-only programme shock,",
         r"on every labour closure. Panels are the labour closures of",
         r"Table~\ref{tab:results-matrix}, rows the utilisation elasticity $\delta$ of",
         r"equation~\eqref{eq:capacity-door}, carried by the goods sectors (CPA A--F), and",
         r"columns the financing closures. Each block reports the change in the price level",
         r"($\Delta P$, the GDP deflator; signed), the largest sectoral price deviation",
         r"($\max|\Delta p|$), employment ($\Delta L$) and household consumption ($\Delta W$),",
         r"all in percent against the reference continuation. The shaded $\delta = 0$ rows",
         r"(no door) reproduce Table~\ref{tab:results-matrix} (the elastic-supply panel at uniform",
         r"$\eta_{s,i} = 0.5$); the door creates no rent. The all-sector profile is in",
         r"Appendix~\ref{app:variants}. Generated from the",
         r"\texttt{pricing\_capacity\_rows-v1-*}, \texttt{pricing\_capacity-v2-*} and",
         r"\texttt{pricing\_mobile-v2-*} run manifests (OpenCGE v0.2.1); no number is",
         r"hand-typed."],
        r"\label{tab:pricing-capacity}",
        [r"% =====================================================================",
         r"% Table 6 --- the capacity door by labour closure and financing (GENERATED)",
         r"% ---------------------------------------------------------------------",
         r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py",
         r"% Goods-only profile (the main-text calibration), all five labour rows",
         r"% (ADR-0032). The all-sector counterfactual lives in tableE_pricing-profiles.tex.",
         r"% ====================================================================="])


def build_mk_table():
    return build_table(
        MK_PANELS, COLS_MK,
        [r"The markup door under the demand-only programme shock, on every labour closure.",
         r"Panels are the labour closures of Table~\ref{tab:results-matrix}, rows the markup",
         r"slope $\kappa$ of equation~\eqref{eq:markup-door} and the routing of the rent",
         r"(Table~\ref{tab:routing}): retained by firms, or distributed to the household,",
         r"which saves the share $s_\pi$ of it ($s$ is the household saving rate). Columns and",
         r"conventions as in Table~\ref{tab:pricing-capacity}; $\Xi/Y$ is the rent's share of",
         r"factor income in percent. $\kappa = 0$ reproduces Table~\ref{tab:results-matrix};",
         r"$\kappa = 0.25$ with distributed rent is in Appendix~\ref{app:variants}. Generated",
         r"from the \texttt{pricing\_\{markup,profit\}-v2-*} run manifests (OpenCGE",
         r"v0.2.1); no number is hand-typed."],
        r"\label{tab:pricing-markup}",
        [r"% =====================================================================",
         r"% Table 7 --- the markup door by labour closure and financing (GENERATED)",
         r"% ---------------------------------------------------------------------",
         r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py",
         r"% Retained (Stage A) and distributed (Stage B) rent, kappa in {0.25, 0.5},",
         r"% all five labour rows (ADR-0032). kappa = 0.25 distributed: tableE.",
         r"% ====================================================================="],
        fit=True)


def build_app_table():
    return build_table(
        APP_PANELS, COLS_MK,
        [r"The pricing-door rows the main text does not carry, on every labour closure.",
         r"Panel A extends the capacity door's $\delta$ to every sector (the all-sector",
         r"counterfactual of Table~\ref{tab:pricing-capacity}); panel B distributes the",
         r"rent of the $\kappa = 0.25$ markup (Table~\ref{tab:pricing-markup} carries",
         r"$\kappa = 0.5$). Conventions as in Table~\ref{tab:pricing-markup}. Generated from",
         r"the same run manifests; no number is hand-typed."],
        r"\label{tab:pricing-profiles}",
        [r"% =====================================================================",
         r"% Table E --- pricing-door profiles not carried in the main text (GENERATED)",
         r"% ---------------------------------------------------------------------",
         r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py",
         r"% Panels: A = capacity door, all-sector counterfactual; B = markup door,",
         r"% kappa = 0.25, rent distributed; each on all five labour rows.",
         r"% ====================================================================="],
        fit=True)


if __name__ == "__main__":
    for path, build in ((OUT_CAP, build_cap_table), (OUT_MK, build_mk_table),
                        (OUT_TABLE_APP, build_app_table)):
        with open(path, "w") as fh:
            fh.write(build())
        print("wrote", path)
