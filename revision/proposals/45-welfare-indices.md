---
title: "Welfare indices for the closure matrix: decisions, specification and workplan"
author: Hermes Agent
date: 2026-10-16
last-updated: October 2026
tags: [welfare, indices, closure-matrix, beyondhulten, revision]
project: BeyondHulten
---

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
- **Section 4.3 wording.** \textcolor{revisionV1}{The sentence claiming the exact measure "corrects a pure quantity index ... for which the former is only a second-order approximation" misdescribes both margins: on prices the approximation is third order under this numeraire, so there is nothing to correct, and under F1 the two differ for a compositional reason, not an approximation one. Suggested replacement: \emph{where prices are pinned the two coincide; where they move the quantity index tracks the exact measure to higher order; and they diverge only under the preference reallocation of F1, where the two are not measuring the same object.}} The paragraph also still carries the typos `assessemnt`, `functioj`, `bthe`, `approaces`, `variabel`, `alredy`, `In what follow`, `gives into estimate`.

## Revision Log

- **Version 1** (October 2026)
- **Version 2** \textcolor{revisionV1}{(October 2026)} — Added the section on why the Tornqvist and the exact measure agree, with the measured order of the residual; sharpened canary C2 and the Section 4.3 wording item.
