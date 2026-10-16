---
title: "Welfare indices for the closure matrix: decisions, specification and workplan"
author: Hermes Agent
date: 2026-10-16
last-updated: October 2026
tags: [welfare, indices, closure-matrix, beyondhulten, revision]
project: BeyondHulten
---

**Version 6** \textcolor{revisionV5}{(October 2026)}
**Version 5** \textcolor{revisionV4}{(October 2026)}
**Version 4** \textcolor{revisionV3}{(October 2026)}
**Version 3** \textcolor{revisionV2}{(October 2026)}
**Version 2** \textcolor{revisionV1}{(October 2026)}
**Version 1** (October 2026)

This note fixes the specification of the three welfare indices announced in
Section 4.3 and to be reported in Section 5, and records the implementation
path. Every number lives in `welfare-panel-v10.md`, which `welfare-panel.py`
regenerates from the run manifests; nothing here is hand-typed.

# Decisions

- **D1 bundle.** The second and third indices are defined on the bundle $c+g$, where $c$ is household consumption and $g$ is public absorption including the programme. Investment and exports are outside: they are exogenous constants in this kernel and carry no welfare content.
- **D2 public valuation.** The weight $\lambda$ on public relative to private absorption is set to 1 throughout. The stance is stated in a footnote in Section 4.3 (inserted 2026-10-16). $\lambda = 0$ returns the pure consumption measure and acts as a check.
- **D3 functional form.** $\sigma = 0.9$, the household's own elasticity, so the aggregation is duality-consistent with the demand system. $\sigma = 1$ (geometric) is the natural comparator to the Tornqvist index.
- **D4 saving.** Both conventions are reported: saving $sE$ as a leak, and as deferred consumption.
- **D5 leisure curvature.** The curvature of the labour--leisure preference is pinned to the closure's own labour-supply elasticity $\eta_s$; the band runs over the curvature parameter.
- **D6 F1.** Under F1 the preference weights $d_i$ are shifted, so an F1 number describes bundle composition rather than a utility comparison across unchanged preferences. F1 stays in the panel, with a footnote.

# The three indices

For any cell, let $E = w\sum_i L_i - T_G(p) - T(p) + F$ be household resources and $E^h = (1-s)E$ the consumption block; $g$ is public absorption at cell prices; $\Pi(p)$ is the consumption price index, which is the numeraire.

- **(a) Tornqvist over consumption.** The quantity index currently reported, $\prod_i (c_i/c_{i0})^{(s_i+s_{i0})/2}$ with expenditure shares $s_i$. It is a superlative approximation to (a').
- **(a') exact money-metric over consumption.** $\Delta E^h / E^h_0$. Exact to all orders because $\Pi(p)=1$, so it needs no deflator in any cell.
- **(b) final-demand index over $c+g$.** The CES aggregate $\left[\sum_i \omega_i^{f} (f_i/f_{i0})^{(\sigma-1)/\sigma}\right]^{\sigma/(\sigma-1)}$ with $f_i = c_i + g_i$ and weights $\omega_i^{f} = f_{i0}/\sum_j f_{j0}$. Equivalent to the value ratio divided by a bundle-specific price index.
- **(c) leisure-augmented money-metric on the $c+g$ base.** To first order, $\Delta W_\ell = \left[\Delta V_f - w\,\Delta N\right]/V_{f0}$, with $V_f = E^h + g$ under the leak convention and $V_f = E + g$ under deferred saving. $w$ is the real wage and $\Delta N$ the employment change. The first-order form is parameter-free: the baseline calibration makes the marginal disutility of work equal to the wage, so the extra wage income and the leisure given up coincide to first order and only the curvature survives.

# Why the Tornqvist and the exact measure agree

\textcolor{revisionV1}{The Tornqvist quantity index is exact for a translog aggregator, not for CES, so at $\sigma = 0.9$ it is \emph{not} the exact money-metric in general. At $\sigma = 1$ (Cobb--Douglas) the two coincide exactly and nowhere else. In this model they nevertheless agree to within $10^{-6}$ on the price-moving cells, and the reason is structural: pinning the CES price index as the numeraire forces the share-weighted log price change to vanish to second order, so the leading Tornqvist--CES discrepancy terms cancel and only a higher-order residual survives. A synthetic check on the model's own structure fits that residual to a cubic in the price dispersion (fitted log-log exponent $3.09$), giving $1.6\times10^{-4}$ pp at the model's dispersion of $0.3$. What makes the term visible is $|\sigma - 1|$, not the dispersion: at $\sigma = 1.5$ and dispersion $0.5$ it reaches $0.02$ pp. The measured null on the F2 and F3 cells is therefore expected, and it checks the reconstruction rather than posing a puzzle.}

\textcolor{revisionV1}{The F1 cells are the exception, and they are not an approximation effect. Their gap is $0.043$ pp, two orders of magnitude above the price residual, and it appears where prices are pinned --- where the approximation gap is exactly zero. It is a compositional artifact of the preference shifter $d_i$: the Tornqvist weights are geared to price-driven substitution, so a pure reallocation of demand weights is not captured by them. The first-order term cancels exactly, because the tilt reallocates expenditure shares without changing their sum, so the leading surviving term is second order in the programme size and is scaled by the programme's concentration relative to the consumption shares. That its size is identical in three cells with different labour closures but the same tilt is the proof of the mechanism.}

\textcolor{revisionV1}{The practical consequence is that the choice between the two indices carries almost no numerical weight --- except under F1, where the two are not answering the same question at all, since the tilt changes the preference ordering.}

# Correction: the final-demand index is a CES aggregate, not the value ratio

\textcolor{revisionV2}{The \texttt{(b)} column of \texttt{welfare-panel-v10.md} is wrong, and the methodological claim it rested on --- ``where prices are pinned, every quantity index reduces to the value ratio'' --- is false. The rule is narrower than that: fixed-weight \emph{arithmetic} indices (Laspeyres, Paasche) do equal the value ratio at unchanged prices, but \emph{geometric} and \emph{curvature-weighted} indices do not, unless every good in the basket grows by the same percentage. The $c+g$ basket does the opposite: private demand scales roughly uniformly while the programme moves seven sectors only, so it is lopsided by construction and the two counts diverge at first order in the reallocation.}

\textcolor{revisionV2}{Recomputed on the actual calibration basket (\texttt{welfare-fd-ces.py}, 62 sectors, private $0.6917$ + public $0.2141$, programme $1.469$\% of the bundle, concentrated in sectors holding $11.2$\% of it):}

| Cell | value ratio (old \texttt{(b)}) | CES $\sigma=0.9$ | CES $\sigma=1$ | CES $\sigma=2$ | Tornqvist | sign flip |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `ALPHA-F1` | +0.0762 | $\mathbf{-0.0339}$ | -0.0238 | +0.0242 | +0.0805 | yes |
| `ALPHA-F2` | +0.0776 | $\mathbf{-0.1425}$ | -0.1229 | -0.0284 | +0.0895 | yes |
| `ALPHA-F3` | +0.0776 | $\mathbf{-0.1425}$ | -0.1229 | -0.0284 | +0.0895 | yes |
| `GAMMA-F1` | $-$0.0946 | $-$0.2042 | $-$0.1941 | $-$0.1464 | $-$0.0904 | no |
| `GAMMA-F2` | +0.2986 | +0.0783 | +0.0979 | +0.1925 | +0.3104 | no |
| `GAMMA-F3` | +3.1986 | +2.9715 | +2.9918 | +3.0894 | +3.2105 | no |

\textcolor{revisionV2}{Three cells change sign. The reported $\sigma=0.5$ figure for ALPHA-F2 reproduces at $-$0.283\%, and at the decided $\sigma=0.9$ the cell reads $-$0.142\% instead of $+0.078\%$. The flip is robust across $\sigma \in [0.5, 2]$ for ALPHA-F2; what $\sigma$ sets is the magnitude, and $\sigma < 1$ (complements) is precisely the regime that penalises a lopsided basket. Note that even the Cobb--Douglas case ($\sigma = 1$, $-$0.123\%) flips: the divergence is not about $\sigma$ alone but about geometric versus arithmetic aggregation.}

\textcolor{revisionV2}{What falls with the claim, beyond the column itself:}

- \textcolor{revisionV2}{The statement that $\sigma$ is irrelevant where prices are pinned. It is exactly backwards: at pinned prices $\sigma$ is the entire content of index (b). Decision D3 needs restating accordingly.}
- \textcolor{revisionV2}{The claim that (b) is computable from the manifests for the nine pinned cells. It is computable for no cell from the manifests; the basket is needed. It can be reconstructed from the calibration artifacts, provisionally, which is what \texttt{welfare-fd-ces.py} does.}
- \textcolor{revisionV2}{The deflator reasoning. No deflator is involved at all: a CES aggregate is computed from the basket directly. A bundle price index would only be needed for the separate deflated-expenditure index, which is a different object, not an equivalent of the CES aggregate.}
- \textcolor{revisionV2}{Table 2 (the $\lambda$ band) and the (c) columns on the $c+g$ base, which were both built as value-ratio families. Table 2 understates the public-valuation effect and must be recomputed as a CES family; the (c) columns on the $c+g$ base are \emph{withdrawn pending re-specification}, because a leisure correction applied to an index expressed in different units is not defined until the aggregation and the numeraire of the correction are stated together.}

\textcolor{revisionV2}{What survives: the consumption-index comparison. Where prices are pinned and the weights are not tilted, consumption quantities scale proportionally, so the Tornqvist, the exact money-metric and the value ratio all agree --- that is why the F2/F3 cells show a null, and it is proportional scaling, not price-pinning, that delivers it. The third-order residual result stands for a CES-demanded bundle and says nothing about the $c+g$ basket.}

# Does the deflator need adapting?

Yes --- a new one, and only for six cells.

- **(b) and (c) need a bundle-specific price index** $\Pi_{c+g}(p)$, built from the $c+g$ baseline shares.
- **The existing seven-component GDP deflator cannot be reused.** It aggregates macro components (including investment, exports, the import and tax leaks) rather than the $c+g$ bundle across sectors.
- **The consumption index $\Pi(p)$ cannot be reused either:** same construction, wrong weights. It is additionally the numeraire, so it is identically 1 and carries no information.
- **Where prices are pinned --- nine of the fifteen cells --- no deflator is needed at all,** because every quantity index then reduces to the value ratio. In those cells the result is also $\sigma$-free, so D3 only binds in the six price-moving rows: the no-reallocation cells and the sectoral BETA column.

# Structural results

At a wage-pinned cell ($w=1$, pinned real wage, prices at baseline) with $\lambda=1$ and deferred saving, the leisure-augmented measure collapses to an exact accounting statement, verified to machine precision on the v10 manifests:

- **F2:** $\Delta V_f - w\,\Delta N = 0$ exactly ($+0.000{\times}10^{0}$ at GAMMA-F2). A tax-financed programme leaves the valued bundle unchanged --- the household pays exactly what the programme costs --- and the extra work is exactly compensated at the margin.
- **F3:** $\Delta V_f - w\,\Delta N = \Delta B_{gov}$ exactly ($+0.013310$ at GAMMA-F3). The externally financed programme delivers exactly the programme's value, as a pure transfer.

These two are the sharpest available canaries; they also state the paper's financing result in one line each.

# Canaries

- **C1 baseline.** Every index equals 1 at the reference continuation.
- **C2 coincidence.** Where prices are pinned and $d_i=1$, the exact money-metric and the Tornqvist index agree to solver precision. Verified at manifest level: $-1.8226$ against $-1.8226$ in the F2/F3 mobile rows. \textcolor{revisionV1}{The agreement is not a knife-edge: it follows from the numeraire and holds to third order in the price change even where prices move.}
- **C3 F1 gap.** The gap in the F1 cells reproduces at a constant $+0.0433$ pp. Regression guard.
- **C4 sigma limit.** As $\sigma \to 1$ the CES aggregate converges to the geometric index.
- **C5 envelope identity.** With $\Delta T_G = \Delta T = \Delta F = 0$ and $s \to 0$, the leisure correction vanishes for any $\Delta N$.
- **C6 time-endowment invariance.** The first-order leisure measure is unchanged when the assumed time endowment is scaled.
- **C7 F2 neutrality.** $\Delta V_f - w\Delta N = 0$ at a wage-pinned F2 cell.
- **C8 F3 transfer.** $\Delta V_f - w\Delta N = \Delta B_{gov}$ at a wage-pinned F3 cell.
- **C9 lambda boundary.** $\lambda = 0$ reproduces index (a') exactly.

# Finding: the bundle indices are not identified under F1

\textcolor{revisionV4}{A user question --- why $W_a$ and $W_b$ disagree under F1, when the programme is already inside household consumption --- turned up a third defect, this time in the F1 rows of the ALT table. Reproduce with \texttt{experiments/welfare-f1-decomp.jl}.}

- \textcolor{revisionV4}{\textbf{The programme really is inside consumption.} \texttt{additive\_demand} returns zeros for any financing other than F2/F3, so under F1 $\sum_i g_i$ is unchanged to machine precision ($\Delta=0$), and the $\lambda$ ladder is pure dilution by the bundle share $C/(C+G)=0.7636$ --- predicted and matched to $0.002$ pp in all five cells.}
- \textcolor{revisionV4}{\textbf{The disagreement is therefore the functional form on the SAME basket.} The CES aggregate over the consumption vector alone ($\lambda=0$) sits $\approx 0.47$ pp \emph{below} the exact money-metric in every F1 cell: $-0.475$, $-0.470$, $-0.474$, $-0.468$, $-0.468$ against $+0.142$, $+0.100$, $+0.420$, $-0.124$, $-0.124$. The regularity across five cells with very different employment and price responses is the diagnosis.}
- \textcolor{revisionV4}{\textbf{Why.} A CES quantity aggregate is the exact utility index only when the observed demands come from its own weights and elasticity. The F1 tilt replaces those weights, so the index loses its foundation and --- with $\sigma<1$ (complements), which is the decided value --- scores the intended reallocation as a loss. It is measuring the policy's purpose as a cost.}
- \textcolor{revisionV4}{\textbf{The contrast is exact.} In F2/F3 the tilt is absent, the private part scales proportionally, and the CES over consumption equals the exact money-metric to the printed digits (gap $0.0000$; e.g. ALPHA-F2 at $-1.823$ on all three). So only in F1 is the bundle index uninterpretable; in F2/F3 the sole difference between $W_a'$ and $W_b$ is the public addition, which is the intended effect.}

\textcolor{revisionV4}{Actions taken (2026-10-16), all reversible:}

- \textcolor{revisionV4}{$W_b$ and $W_c$ are suppressed on the F1 rows of the ALT table, marked $\ddag$ at the group header and explained in its caption.}
- \textcolor{revisionV4}{$W_a'$ is the F1 head-line, carrying the composition caveat rather than a welfare reading.}

\textcolor{revisionV4}{\textbf{Does a preference-consistent re-weighting rescue the F1 cells?} A follow-up in the same round, and a correction to the paragraph above: \emph{it depends on the convention}, and the two coherent readings do not agree. Read as \emph{each state under its own preferences}, in money-metric units (weights normalised to sum to one, so $U=E$ at unit prices), the bundle index collapses \emph{identically} onto the real-income measure, because a CES household's utility \emph{is} its real income: on the fixture $U_\omega(c_0)=1.000000000000$ and $U_\alpha(c_1)=1.050000000000$, ratio $1.050000000000$ --- the income effect to twelve digits. Under that reading $W_b$ would indeed equal $W_a'$, and the column would be a duplicate. Read instead as \emph{one updated aggregator applied to both} bundles --- is the old basket any good by the new tastes? --- it does \emph{not} collapse: $+24.6$~\% against a $+5.0$~\% income effect, with the baseline-weighted aggregate at $-16.0$~\%, a spread of some 40~pp around 5~pp. \texttt{welfare-preference-shift.py} reproduces all three. The one object that is consistent, free of cardinalisation and about the household is $E^h/\Pi(p)$ --- and $\Pi(p)=1$ on every solve path, verified to $2.8\times10^{-12}$ across all fifteen cells (\texttt{experiments/welfare-cpi-check.jl}) --- which is exactly $W_a'$. So the F1 bundle column is either a duplicate of $W_a'$ or a cross-regime ranking that is not a statement about the household's welfare. Neither earns a column.}

\textcolor{revisionV4}{\textbf{The identity behind the collapse, and the three candidate numbers.} Recovering the preference weights from the solutions themselves ($\alpha_i \propto c_i p_i^{\sigma}$, exact at unit prices) shows why the `each state by its own weights' reading collapses: with own weights the CES aggregator reduces identically to $\sum_i q_i$, so evaluating a state under its own weights is nothing but \emph{that object's value ratio}. Applied to consumption it is therefore exactly $W_a'$ --- confirmed on the real 62-sector cells, $+0.0999$ against $+0.0999$ (ALPHA-F1) and $-0.1238$ against $-0.1238$ (GAMMA/DELTA-F1). Applied instead to the bundle $c+g$ it is the value ratio $+0.076$, differing from $W_a'$ purely by the $23.6$~\% unchanged public share. So an F1 row has three candidate numbers for $W_b$, each defensible under a stated convention: $+0.076$ (own weights $\equiv$ value ratio), $-0.281$ (baseline weights), $+0.690$ (one updated aggregator applied to both bundles). The first equals $W_a'$ only if $g$ is dropped from the argument --- i.e. only by no longer being $W_b$. \texttt{experiments/welfare-f1-reweight.jl} reproduces all three. The gap between the per-state reading and the single-updated-aggregator reading is definitional, not empirical: the two share a numerator --- the counterfactual bundle valued by the counterfactual preferences, $0.69236601$ --- and differ only in the \emph{denominator}, the valuation of the \emph{baseline} basket: $0.69167530$ by its own weights against $0.68761994$ by the updated ones, a revaluation of the benchmark by $-0.59$~\% (ALPHA-F1). The updated aggregator therefore inflates the F1 result by re-pricing the status quo, not by measuring anything about the counterfactual. Note the incidental result: under F1 the withdrawn value ratio \emph{is} the preference-consistent index for the bundle, while under F2/F3 it is the CES --- the two regimes genuinely admit different consistent conventions, which is the whole reason the F1 column cannot be made commensurable with the others.}

\textcolor{revisionV4}{One related precision: F1 does \emph{not} imply $\eta_s=0$. F1 is a financing closure and the labour closure still sets the supply elasticity, so only the BF and ALPHA rows have an inert leisure term under F1; BETA ($\eta_s=0.5$) and GAMMA/DELTA move employment there ($\Delta L=+0.12,-0.10,-0.10$) and $W_c$ retains a live leisure term.}
- \textcolor{revisionV4}{Table 5's final-demand row is scoped to the F2/F3 programme cells only.}
- \textcolor{revisionV4}{Section 4.3 carries the scope caveat once, at the point where the three indices are announced, rather than cell by cell.}

## Specification: the own-weights bundle index \textcolor{revisionV5}{\normalsize [added v6]}

\textcolor{revisionV5}{\textbf{Decision (operator, 2026-10-16).} The bundle column carries the \emph{own-weights} aggregate rather than the fixed-weight CES: it is the only bundle measure that is preference-consistent in every regime, and it is the number that belongs in the F1 rows.}

\textcolor{revisionV5}{\textbf{Definition.} For a state $s$ with bundle $f_s = c_s + \lambda g_s$ ($\lambda = 1$), let $w_s = f_s/\sum_i f_{si}$ be that state's own bundle-composition weights. The index is $B^{\mathrm{own}} = U(f_1, w_1)/U(f_0, w_0)$ with $U(q,w) = [\sum_i w_i^{1/\sigma} q_i^{(\sigma - 1)/\sigma}]^{\sigma/(\sigma - 1)}$ and $\sigma = 0.9$. \textbf{Property:} with own weights the aggregator reduces identically to $\sum_i q_i$, so $B^{\mathrm{own}} = \sum_i f_{1i} / \sum_i f_{0i}$, the ratio of bundle quantity sums. At base prices this coincides exactly with the bundle's value ratio; where prices move it does not.}

\textcolor{revisionV5}{\textbf{Values, all fifteen cells} (percent against the reference continuation):}

| financing | labour | own-weights | value ratio | fixed-weight CES | $W_c^{leak}$ | $W_c^{def}$ |
| --- | --- | --- | --- | --- | --- | --- |
| F1 | BF | -0.013 | +0.030 | -0.253 | -0.013 | -0.013 |
| F1 | ALPHA | **+0.076** | +0.076 | -0.281 | +0.076 | +0.078 |
| F1 | BETA | +0.261 | +0.293 | -0.040 | +0.129 | +0.149 |
| F1 | GAMMA | -0.095 | -0.095 | -0.451 | +0.013 | -0.000 |
| F1 | DELTA | -0.095 | -0.095 | -0.451 | +0.013 | -0.000 |
| F2 | BF | -0.028 | +0.039 | -0.318 | -0.028 | -0.210 |
| F2 | ALPHA | +0.078 | +0.078 | -0.291 | +0.078 | -0.101 |
| F2 | BETA | +0.272 | +0.310 | -0.059 | +0.133 | -0.027 |
| F2 | GAMMA | +0.299 | +0.299 | -0.069 | +0.160 | 0.000 |
| F2 | DELTA | +0.299 | +0.299 | -0.069 | +0.160 | 0.000 |
| F3 | BF | -0.028 | +0.039 | -0.318 | -0.028 | -0.210 |
| F3 | ALPHA | +0.078 | +0.078 | -0.291 | +0.078 | -0.101 |
| F3 | BETA | +0.272 | +0.310 | -0.059 | +0.133 | -0.027 |
| F3 | GAMMA | +3.199 | +3.199 | +2.827 | +1.234 | +1.331 |
| F3 | DELTA | +3.199 | +3.199 | +2.827 | +1.234 | +1.331 |

\textcolor{revisionV5}{\textbf{The leisure-augmented measure rides the same bundle.} $W_c$ is the first-order leisure-augmented money-metric on the \emph{own-weights} bundle, $W_c = (B^{\mathrm{own}} - 1) - w\,\Delta N / V_0$, so the bundle column and $W_c$ share one parent and one regime-consistency. This is a substantive change at the wage-pinned pole: GAMMA/DELTA under F1 moves from $-0.344$ (CES parent) to $+0.013$ for the leak convention, because employment \emph{falls} there and the leisure term is a gain. BF and the two price-moving rows are treated identically, since the own-weights aggregate is defined without reference to the price level.\footnote{The kernel exposes this as \texttt{welfare\_leisure(...; aggregator = :ownweights)}; the fixed-weight CES variant remains available as \texttt{:ces} and stays the default, so the stanced reading is still computable.}}

\textcolor{revisionV5}{\textbf{Three specification points.} (i) In every price-pinned cell the own-weights aggregate and the bundle's value ratio agree to four decimals; the two price-moving rows (BF, BETA) diverge, and the table reports the own-weights aggregate with the value ratio shown beside it. (ii) This measure \emph{is} the value ratio of the withdrawn index (b) in all price-pinned cells, so it reinstates that number as the bundle column and supersedes the fixed-weight CES there. The CES keeps its separate role as the \emph{stanced} measure --- it values public absorption through private preferences --- and remains the right instrument where the basket moves exogenously rather than by preference. (iii) $W_c$ inherits the same parent (below), so every bundle-based column is regime-consistent.}

\textcolor{revisionV5}{Reproduce with \texttt{experiments/welfare-ownweights.jl} and, for the leisure columns, \texttt{experiments/make\_table4\_alt.jl}.}

# Phase 1 --- implemented

\textcolor{revisionV3}{Landed in \texttt{opencge} (branch \texttt{pricing-door}, uncommitted):}

- \textcolor{revisionV3}{\texttt{src/core/diagnostics.jl}: \texttt{final\_demand\_bundle(model, p, y, w; lambda, external\_transfer)} returning the per-sector private, public and valued bundles; \texttt{ces\_quantity\_index(q, q0, sigma)} with the Cobb--Douglas limit and the zero-baseline guard; \texttt{real\_consumption\_ev(sol, base)}; \texttt{real\_final\_demand\_ces(sol, base; lambda, sigma)}; \texttt{welfare\_leisure(sol, base; lambda, sigma, saving)} returning the index, the value ratio, the employment change, the leisure cost and the measure. All four exported.}
- \textcolor{revisionV3}{\texttt{tests/test\_welfare\_indices.jl}, registered in \texttt{tests/runtests.jl}.}

\textcolor{revisionV3}{Canaries carried: C10 lopsided basket (bread/cement: CES at $\sigma=0.9$ reads $1.048520$ against a value ratio of $1.136364$; every $\sigma<1$ lands below it, and the T\"ornqvist lands above at $1.149083$ --- so the divergence is not repairable by swapping the functional form); proportional growth makes every index agree with the value ratio at every $\sigma$; $\sigma=1$ is the geometric index; the zero-baseline component is dropped rather than divided by zero; the baseline normalisation $Q=1$; $\lambda=0$ returns the consumption measure; the error paths.}

\textcolor{revisionV3}{Evidence: the full suite runs green --- 490 assertions, 0 failures, 0 errors. The C10 testset is the regression guard for the mistake the section above records.}

\textcolor{revisionV3}{Not done: the leisure columns remain first-order only; the bundle deflator is unnecessary (the CES aggregate is computed from the basket, see the correction above), so that item is closed as withdrawn; and the real-data indices still come from the provisional calibration-side reconstruction rather than the kernel.}

# Workplan

- **Phase 0 --- decisions.** Complete; recorded above.
- **Phase 1 --- kernel.** Add `final_demand_bundle`, `real_final_demand_ces(sol, base; bundle, sigma)` and `welfare_leisure(sol, base)`, plus the bundle deflator $\Pi_{c+g}(p)$, in `opencge/src/core/diagnostics.jl`. Keep `real_consumption` as the Tornqvist. Ship canaries C1--C9 in the same commit. Acceptance: `julia --project=. tests/runtests.jl` green, C1--C9 passing.
- **Phase 2 --- harness.** Add `welfare_tq`, `welfare_ev`, `welfare_fd_ces`, `welfare_fd_ces_defl`, `welfare_leisure_leak`, `welfare_leisure_deferred` to `[metrics]`; bump `schema_version` to 3; update `experiments/README.md` and `tests/test_run_manifest.jl`.
- **Phase 3 --- generation.** Mint `matrix_5x3_v11`: the 33 v10 cells plus the uniform $\eta_s = 1.5$ cells, with the new metrics and restored `solution.csv`. Preregister the design hash; v10 remains history (ADR-0004).
- **Phase 4 --- tables.** Keep Table 4 as it stands and add a welfare panel (rows = the fifteen cells, columns = the indices), plus the $\eta_s$ ladder table. Both generated, no hand-typed numbers.
- **Phase 5 --- text.** Section 4.3 states each index and the stance it embeds; Section 5 reports the three side by side and uses C2 to explain why the choice is immaterial where prices are pinned.

# Open items

- **(b) and (c) on the $c+g$ base are blocked for six cells** until $\Pi_{c+g}(p)$ exists in the kernel. Those cells are currently `n/a` in the panel.
- **BF leisure valuation** uses the wage-bill-weighted aggregate wage; it should use the binding sectoral wage.
- **Second-order leisure term** needs the curvature mapping $\eta_s \to \theta$ and the time-endowment convention; the first-order measure is free of both.
- **C10 lopsided basket.** \textcolor{revisionV2}{A CES aggregate must differ from the value ratio whenever the basket does not scale proportionally. Regression guard: ALPHA-F2 at $\sigma=0.9$ reads $-0.1425$, not $+0.0776$. Include a proportional-growth fixture where the two agree exactly.}
- **C7 and C8 scope.** \textcolor{revisionV2}{These two identities hold for the \emph{value} of the bundle. Whether an analogue holds for the CES aggregate is open; do not present them as canaries for index (b).}
- **Section 4.3 wording.** \textcolor{revisionV1}{The sentence claiming the exact measure "corrects a pure quantity index ... for which the former is only a second-order approximation" misdescribes both margins: on prices the approximation is third order under this numeraire, so there is nothing to correct, and under F1 the two differ for a compositional reason, not an approximation one. Suggested replacement: \emph{where prices are pinned the two coincide; where they move the quantity index tracks the exact measure to higher order; and they diverge only under the preference reallocation of F1, where the two are not measuring the same object.}} The paragraph also still carries the typos `assessemnt`, `functioj`, `bthe`, `approaces`, `variabel`, `alredy`, `In what follow`, `gives into estimate`.

## Revision Log

- **Version 1** (October 2026)
- **Version 2** \textcolor{revisionV1}{(October 2026)} — Added the section on why the Tornqvist and the exact measure agree, with the measured order of the residual; sharpened canary C2 and the Section 4.3 wording item.
- **Version 3** \textcolor{revisionV2}{(October 2026)} — Correction: the final-demand index is a CES aggregate, not the value ratio. Recomputed column, three sign flips, four dependent claims withdrawn; canaries C7/C8 rescoped and C10 added.
- **Version 4** \textcolor{revisionV3}{(October 2026)} — Phase 1 implemented in the kernel: five functions plus the canary testset; full suite green at 490 assertions.
- **Version 5** \textcolor{revisionV4}{(October 2026)} — F1 finding: the bundle indices are not identified under a preference change (a ~0.47 pp aggregation artefact on the same basket). F1 rows of the ALT table suppressed; Table 5 scoped; Section 4.3 caveat added.
- **Version 6** \textcolor{revisionV5}{(October 2026)} — Adopted the own-weights bundle index as the bundle column, carried in all rows, and rebased the leisure-augmented measure $W_c$ on the same bundle (kernel gains \texttt{welfare\_leisure(...; aggregator)} plus canary C11, suite 490 $\rightarrow$ 496 assertions); specification, the identity and the fifteen values added as a new section.
