#!/usr/bin/env python3
"""Regenerate the pricing-door tables from the Phase 7 run manifests
(ADR-0029, ADR-0030, ADR-0031, ADR-0032).

No number is hand-typed: every value is a manifest field, except Table 6's
W_D, which no manifest records; it is read from table6_wd.csv, written by
experiments/make_table6_wd.jl from manifest-gated re-solves (as Table 5's W_D).
Run from anywhere:
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
import csv, tomllib, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WD_CSV = os.path.join(REPO, "revision/tables/table6_wd.csv")
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


def load_wd():
    """W_D per run id (percent), from the manifest-gated sidecar of make_table6_wd.jl."""
    with open(WD_CSV, newline="") as fh:
        return {r["run_id"]: float(r["w_d"]) for r in csv.DictReader(fh)}


WD = load_wd()


def values(rid):
    m = manifest(rid)["metrics"]
    return dict(P=100 * (m["gdp_deflator"] - 1.0), D=100 * m["max_abs_price_dev"],
                L=100 * (m["employment"] - 1.0), W=100 * m["consumption_rel"],
                WD=WD.get(rid, float("nan")), S=100 * m.get("profit_share", 0.0))


def num(x, signed=True):
    assert x == x, "missing value (rerun experiments/make_table6_wd.jl?)"
    if abs(x) < 5e-3:
        return r"$0.00$"
    if not signed:
        return "$%.2f$" % x
    return ("$+%.2f$" % x) if x > 0 else ("$%.2f$" % x)


COLS_CAP = [("P", r"$\Delta P$", True), ("D", r"$\max|\Delta p|$", False),
            ("L", r"$\Delta L$", True), ("WD", r"$\Delta W_D$", True)]
COLS_MK = [("P", r"$\Delta P$", True), ("D", r"$\max|\Delta p|$", False),
           ("L", r"$\Delta L$", True), ("W", r"$\Delta W$", True), ("S", r"$\Xi/Y$", False)]


def hdr_cell(h):
    """One column-header cell: centred (the data columns are r-aligned) and, for
    max |dp|, split over two lines so the header stays narrow and each value
    sits under a short label. \\shortstack is base LaTeX (no package needed)."""
    if r"\max|\Delta p|" in h:
        return r"\multicolumn{1}{c}{\shortstack{$\max$\\$|\Delta p|$}}"
    return r"\multicolumn{1}{c}{" + h + "}"


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
    sub = " & ".join(hdr_cell(h) for _, h, _ in cols)
    A(r"\multicolumn{1}{c}{Labour closure and door} & " + " & ".join([sub] * len(FIN)) + r" \\")
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


def build_guided_table(panels, cols, stub, caption_lines, label, comment_lines, fit=False):
    """The pricing-door layout tuned for reading across many columns: financing
    blocks named as in Table 4 and set apart by wider gaps, labour panels separated
    by vertical space instead of rules, and the first row of every panel shaded
    as a guide line across the row (the delta = 0 reference in Table 6 and the
    retained kappa = 0.25 case in Table 7); thin grey rules between the columns
    of a block keep each value under its header (they stop at the panel headings
    and block gaps). The block gaps are empty spacer columns rather than @{}
    glue, and the outer columns keep their padding, so the shading runs unbroken
    from rule end to rule end."""
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
    if fit:
        A(r"\resizebox{\textwidth}{!}{%")
    guide = r"!{\color{tableline}\vrule width 0.3pt}"
    A(r"\begin{tabular}{l" + "".join([r"@{\hspace{3pt}}p{2pt}@{\hspace{3pt}}" + guide.join("r" * k)
                                      for _ in FIN]) + r"}")
    A(r"\toprule")
    A(" & " + " & ".join([r"& \multicolumn{%d}{c}{\textbf{%s}}" % (k, f) for f in FIN]) + r" \\")
    A(" & " + " & ".join([r"& \multicolumn{%d}{c}{\scriptsize %s}" % (k, FIN_NOTE[f]) for f in FIN])
      + r" \\")
    A("".join([r"\cmidrule(lr){%d-%d}" % (3 + (k + 1) * i, 2 + (k + 1) * i + k)
               for i in range(len(FIN))]))
    sub = " & ".join(hdr_cell(h) for _, h, _ in cols)
    A(r"\multicolumn{1}{c}{" + stub + "} & " + " & ".join(["& " + sub] * len(FIN)) + r" \\")
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
    A(r"\end{tabular}" + ("}" if fit else ""))
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
        [r"Integrating a Kaldor--Verdoorn term in the general closure matrix. Each block reports",
         r"the change in the price level ($\Delta P$, the GDP deflator; signed), the largest sectoral",
         r"price deviation ($\max|\Delta p|$), employment ($\Delta L$) and the quantity index of demand",
         r"over $c+g$ ($\Delta W_D$, see also Table~\ref{tab:welfare-comparison}), all in percent against",
         r"the reference continuation. The shaded $\delta = 0$ rows (no door) reproduce",
         r"Tables~\ref{tab:results-matrix} and~\ref{tab:welfare-comparison}."],
        r"\label{tab:pricing-kaldor}",
        [r"% =====================================================================",
         r"% Table 6 --- the capacity door by labour closure and financing (GENERATED)",
         r"% ---------------------------------------------------------------------",
         r"% Do NOT edit by hand. Regenerate with:  make_pricing_door.py",
         r"% Goods-only profile (the main-text calibration), all five labour rows",
         r"% (ADR-0032). The all-sector counterfactual lives in tableE_pricing-profiles.tex.",
         r"% ====================================================================="])


def build_mk_table():
    return build_guided_table(
        MK_PANELS, COLS_MK, "Markup and rent routing",
        [r"The markup door under the demand-only programme shock, on every labour closure.",
         r"Panels are the labour closures of Table~\ref{tab:results-matrix}, rows the markup",
         r"slope $\kappa$ of equation~\eqref{eq:markup-door} and the routing of the rent",
         r"(Table~\ref{tab:routing}): retained by firms, or distributed to the household,",
         r"which saves the share $s_\pi$ of it ($s$ is the household saving rate). Columns and",
         r"conventions as in Table~\ref{tab:pricing-kaldor}, except that $\Delta W$ is household",
         r"consumption ($W_C$ of Table~\ref{tab:results-matrix}); $\Xi/Y$ is the rent's share of",
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
         r"counterfactual of Table~\ref{tab:pricing-kaldor}); panel B distributes the",
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
