---
title: "Variation in GAMMA: options for generating outcome variation in the fixed-wage closure"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [assessment, labour-closure, gamma, wage-structure, options, evidence]
last-updated: September 2026
---

**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

This document is an option catalogue with measurements, not a plan: adopting any
option below requires an ADR and a registry entry, per the house rules. Closure
definitions and their current state live in `registry/closures.toml` and
`docs/status.md`; they are not restated here (ADR-0003).

# The problem: GAMMA carries one degree of freedom

The fixed-wage closure pins one real wage as the numeraire, keeps the allocation
fully mobile and drops both the labour equation and the external unknown. Its
only free dimension in the evaluation matrix is therefore the financing row, and
that row is nearly degenerate:

| Source of variation in GAMMA | Cells that differ | Measured spread |
| --- | --- | --- |
| Financing F1 (compositional tilt) | `GAMMA-F1` | `gdp_rel` -0.000973, `consumption_rel` -0.000806 |
| Financing F2 (tax) and F3 (external) | `GAMMA-F2`, `GAMMA-F3` | identical real allocation (financing neutrality), `gdp_rel` +0.001260, `consumption_rel` -0.015333 |
| Labour closure: DELTA | `DELTA-F1..F3` | reproduces GAMMA to six digits (Leontief corner at `epsilon = 1e-4`) |

So of the five labour rows, GAMMA and DELTA contribute essentially one degree of
freedom, and within GAMMA the F2 and F3 columns are real-neutral by theorem. The
flatness is not a property of the model class: it is a consequence of the design
setting the wage vector to its degenerate value, all sectors equal. The kernel
then hard-codes that value as the scalar numeraire.

# Where variation can come from

Four doors, ordered by cost and by how much of the closure they disturb.

## Door 1: the wage structure (pinned wage vector)

Replace the scalar pin `w = 1` by a pinned vector `wbar`. The system keeps its
2N unknowns and 2N equations; the wage enters as a parameter rather than an
unknown. This is the door developed in the rest of this document: it reuses the
vector-wage apparatus built for the eta = 0 endpoint (ADR-0020), it needs no new
closure, and it is the only door that adds a *continuous* dimension.

## Door 2: sectoral supply shocks

Already an input of the design schema, no code change: a sectorally concentrated
cost shifter moves prices with the wage pinned, so real wages, employment and the
external position all move. This door also makes the allocation margin bite,
which is why the assessment records that `eta_s` is identifiable only under a
supply-side scenario. It is the cheapest way to widen GAMMA's outcome range, and
it is complementary to Door 1 rather than a substitute: Door 1 varies the wage
structure at a given shock, Door 2 varies the shock at a given wage structure.

## Door 3: segmented wage setting

Pin the wages of a subset of sectors and let the rest share one flexible wage:
unknowns `[p; y; w_free]` (2N+1), with the flexible side carrying a labour-market
equation. This nests GAMMA (all sectors pinned) and ALPHA (none pinned) and makes
the sticky share a continuous parameter. It is the most informative variant for a
labour-market story, but it is a new closure with its own admission question, and
it needs its own ADR.

## Door 4: nominal wage rules (not proposed here)

A fixed *nominal* wage (real wage falling with the CPI) is a different closure: it
needs a nominal anchor other than the wage, which the current model does not
have. `docs/ETAs.md` records a no-go on aggregate nominal wage rules; that record
was not re-read for this document and Door 4 is therefore not proposed. If it is
ever pursued, it must be checked against that no-go first, and the indexation
degree would be its parameter.

# The wage-structure door in detail

## What changes in the system

The residual is the existing fixed-wage system with the scalar `1.0` replaced by
`wbar`: `p = cost(p, wbar)` and `y = demand(p, wbar)` for every sector, with the
household wage income `sum(wbar_i * L_i)` and the CPI numeraire implicit in the
wage vector itself. At `wbar = 1` it is the current GAMMA system, and the
prototype reproduces the `GAMMA-F2` cell (residual 8.9e-16, employment
1.00125981) from the same warm start the harness uses.

## What is not an instrument: the wage level

A uniform rescale of the wage vector is a numeraire change, not a shock. The unit
cost block is homogeneous of degree one in `(p, w)`, and the demand block is
homogeneous of degree zero in `(p, w)` because nominal income scales with `w` and
the CPI with `p`; hence `(c * p*, y*)` solves the system for the pin `c * wbar`,
and by determinacy it is the solution. Measured on the `GAMMA-F2` cell:

| Pin | max abs log y/y0 | Employment | CPI | Real income |
| --- | ---: | ---: | ---: | ---: |
| `wbar = 1` | 0.0e+00 | 1.00125981 | 1.00000000 | 1.00125981 |
| `wbar = 0.9` | 6.7e-16 | 1.00125981 | 0.90000000 | 1.00125981 |
| `wbar = 1.1` | 8.9e-16 | 1.00125981 | 1.10000000 | 1.00125981 |

Two consequences. First, a statement about wage *restraint* (a level statement)
is outside this closure's language: the closure can only speak about the wage
*structure*. Second, the invariance is a property of this model class, not a
general truth: it holds because there is no nominal anchor besides the wage. With
a nominal government bill, nominal debt or money, a uniform wage change would be
real.

## What is an instrument: the wage structure

A relative change in the pinned vector is a genuine shock. Measured on the same
cell, all residuals at or below 2e-13:

| Scenario | max abs log y/y0 | Employment | CPI | Real income |
| --- | ---: | ---: | ---: | ---: |
| `wbar = 1` (baseline) | 0.0e+00 | 1.00125981 | 1.00000000 | 1.00125981 |
| sector 1 +10 percent | 2.1e-02 | 1.00116732 | 1.00092606 | 1.00080639 |
| sector 1 -10 percent | 2.3e-02 | 1.00139321 | 0.99902004 | 1.00174815 |
| programme sectors +10 percent | 1.6e-02 | 1.00686452 | 1.00542537 | 1.01186464 |
| programme sectors -10 percent | 1.7e-02 | 0.99589845 | 0.99426945 | 0.99040765 |
| sectors 1 to 5 +10 percent | 3.5e-02 | 1.00116698 | 1.00128008 | 1.00063736 |

A ten percent wage push in the seven programme sectors moves employment by about
half a percent and real income by about one percent, which is larger than
GAMMA's entire existing F1 to F2 spread. The instrument is a continuous
N-dimensional parameter, so the cell count of the matrix stops bounding the
number of distinct GAMMA economies.

## The code surface a promotion touches

- `problem_fixed` and `_solve_fixed` need the wage vector threaded through as a
  parameter (the scalar call sites keep their arithmetic).
- `external_balance_canary` and `gdp_components` hard-code `w = 1` on their 2N
  branch. \textcolor{revisionV1}{Version 2 measured the identity through the price-weighted clearing residual instead (validated against the canary along quantity perturbations), so the closure question is answered; the extension is still required so that the acceptance gate evaluates a tilted pin on the kernel's own arithmetic rather than on a proxy.}
- The harness reports the scalar `wage` metric; as at eta = 0 it should be the
  wage-bill-weighted aggregate, with `wage_min` and `wage_max` carrying the
  structure.
- The F3 booking is priced at the equilibrium prices, so it moves with the wage
  structure; the identity gate should be re-checked per scenario.

## Determinacy and conditioning

The admission criterion of the fixed-wage regime involves only technology and
demand shares (`A_bill`, `lambda`, `m`, `s`) and never `w`, so it is invariant
across the wage-structure scenarios; the measured Jacobian condition number is
flat (6.32 to 6.34) and the residual floor is unchanged. The fixed-wage system is
well conditioned at this calibration, which also settles a wording question from
the v6 work: the GAMMA/DELTA metric sensitivity to the warm start is a
measurement-layer amplification of the residual level, not a singular solve.

# Should these numbers be reproduced?

Yes, and the claims should be split before anything is cited.

| Claim | Strength today | What reproduction adds |
| --- | --- | --- |
| Uniform rescale is a numeraire change | Derived from the block's homogeneity and confirmed numerically to 1e-15 | Nothing material; keep as a property statement |
| A tilted vector moves quantities, employment and prices | Measured, but on a prototype that *mirrors* the fixed-wage residual | The kernel path, because translation is where defects appear |
| The external account closes under a tilted vector | \textcolor{revisionV1}{Measured in v2: the price-weighted clearing residual stays at the solver-residual level (at most 7.0e-12), and it equals the kernel canary along quantity perturbations} | \textcolor{revisionV1}{Nothing further for the closure; the canary must still take the vector so that the gate runs on the kernel's own arithmetic} |
| The magnitude of the responses | Prototype numbers, scenario-dependent | A pre-registered design, so the scenarios are fixed before the run |

\textcolor{revisionV1}{Version 2 closes the decisive gap. The identity under a tilted vector is now measured: the price-weighted clearing residual, which is the quantity the canary computes, stays at the solver-residual level under every scenario (2.9e-16 at the degenerate pin, at most 7.0e-12 under a tilt). That proxy was validated against the kernel's own canary off equilibrium, where the gap is large: the two agree exactly along quantity perturbations (ratio 1.0000 at one and ten per mille), which is why the check is credible even though both quantities are machine zero on the equilibrium itself. One limit remains, recorded in the open items: the relation is an equilibrium statement, not a universal algebraic identity (a price perturbation gives -1.95), so the canary should still be extended before the gate is relied upon for a tilted pin. ADR-0021 (proposed) carries the adoption decision.}

The cautionary precedent is from this same work: the eta = 0 prototype was
correct, and the defect appeared when it was translated into the kernel, where
the first draft used the frozen allocation in the first-order condition, making N
equations identically zero. The external-account gate caught it within one run.
That is the argument for reproducing rather than quoting: a mirrored residual is
not the residual the harness will run, and the gate that protects the paper runs
on the harness.

The reproduction is cheap and follows the established path: promote (dispatch
only, then assert that the twelve non-BF cells and the GAMMA/DELTA cells
reproduce their v6 manifests), extend the canary and `gdp_components`, register a
handful of wage-structure scenarios in a design, preregister it, and run it as a
new generation. Only the last step produces citable numbers.

# Do the numbers rest on an ad-hoc assumption?

No on the device, and yes on the scenario values, which is the normal status of a
shock in this repository. The classification:

| Element | Status | Why |
| --- | --- | --- |
| Pinning a real wage | Not ad-hoc | It is the existing fixed-wage closure; the vector is its natural generalization, not a new device |
| The vector being a parameter | Not ad-hoc, but exogenous | The wage structure is an institutional scenario, like the programme incidence or the shock size; it is not an equilibrium object |
| Level invariance | Closure-conditional | It follows from the absence of a nominal anchor besides the wage; it should be stated with that condition |
| The magnitudes | Scenario-dependent | They are instrument responses, not estimates; report the mapping over pre-registered scenarios, never a single headline number |
| Technology and demand parameters | Inherited | `theta`, `epsilon`, `sigma`, `eta_s` and the calibration are unchanged from the matrix |

The comparison that matters is with the retired allocation wedge, which was
dropped precisely because it was never derived as a coefficient of the model. The
wage vector is the opposite case: the parameter is already in the closure's
language (the fixed wage is its degenerate point), and the current design simply
sets it to that point. What must not happen is the reverse inference: that a
wage-structure scenario identifies a labour-supply elasticity, or that a
welfare difference across scenarios is a statement about wage policy. It does
not, and the document claims no such thing.

# What the paper gains

The two labour margins become separable instead of being conflated in one
endpoint:

| | allocation flexible | allocation frozen |
| --- | --- | --- |
| wage flexible | ALPHA, BETA | (retired interpolated wedge) |
| wage pinned as a vector | GAMMA with a wage structure (this document) | BF at eta = 0 |

GAMMA with a wage structure and BF then differ in exactly one respect, the
allocation margin, which is the comparison the matrix cannot currently make:
GAMMA's wage vector is degenerate and, before ADR-0020, BF was a normalization of
ALPHA. A DELTA companion (the Leontief corner with a wage structure) is a cheap
extension that would test whether the corner amplifies pass-through.

# Options and provenance

Per option, what it needs and where the evidence stands.

| Option | What it needs | Evidence today |
| --- | --- | --- |
| Door 1, wage structure | Thread `wbar` into `problem_fixed` and `_solve_fixed`; extend the canary and `gdp_components`; new closure id, ADR, design | `experiments/probes/probe11_gamma_wage_structure.jl`; the measurements above; ADR-0020 for the vector-wage apparatus |
| Door 2, sectoral supply shocks | A design with sectoral supply shocks; no code change | The shock schema in `experiments/README.md`; the assessment's identifiability note on `eta_s` |
| Door 3, segmented wages | A new closure: pinned subset plus one flexible wage; ADR and admission check | None yet; sketched here only |
| Door 4, nominal wage rules | A nominal anchor the model does not have; check `docs/ETAs.md` no-go first | `docs/ETAs.md` (not re-read here) |

Provenance of the claims in this document, per claim:

| Claim | Source |
| --- | --- |
| F2 and F3 are real-neutral in GAMMA | Financing neutrality, `runs/matrix_5x3-v6-GAMMA-F*` |
| DELTA reproduces GAMMA to six digits | The v6 flow table, `paper/tables/matrix_5x3_v6_flows.md` |
| The baseline pin reproduces the kernel cell | `probe11`, residual 8.9e-16, employment 1.00125981 |
| Level invariance | Homogeneity argument plus the uniform-rescale rows of `probe11` |
| The structure moves outcomes | The tilted-vector rows of `probe11` |
| The identity closes under a tilted vector | \textcolor{revisionV1}{Measured in v2: the price-weighted clearing residual (2.9e-16 degenerate, at most 7.0e-12 tilted), validated against the kernel canary along quantity perturbations} |
| \textcolor{revisionV1}{The canary cannot be used directly for a tilted pin} | \textcolor{revisionV1}{\texttt{external\_balance\_canary} pins w = 1 on its 2N branch; hence the proxy} |
| The conditioning is flat | `probe11`, `cond(J)` in 6.32 to 6.34 |

# Evidence and open items

Evidence: `experiments/probes/probe11_gamma_wage_structure.jl` (this document's
measurements), `probe7_sectoral_wages_eta0.jl` (the eta = 0 closure the apparatus
comes from), `probe8_promotion_verification.jl` and `probe9_nonbf_reproduction.jl`
(the promotion's verification pattern to copy). Runs: `runs/matrix_5x3-v6-GAMMA-*`
and `runs/matrix_5x3-v6-DELTA-*`.

Open items, in the order a promotion would meet them:

- \textcolor{revisionV1}{The identity under a tilted vector is measured (v2), so this item no longer blocks: what remains is that the relation used to measure it is an equilibrium statement, not a universal algebraic identity (a price perturbation gives -1.95), and that the kernel canary pins w = 1 on its 2N branch. A promotion must extend the canary and gdp\_components so that the acceptance gate runs on kernel arithmetic, and must re-assert the identity per scenario.}
- The wage-structure scenarios must be fixed by preregistration before the run,
  otherwise the variation is tuned rather than reported.
- `docs/ETAs.md` carries a no-go on aggregate nominal wage rules that was not
  re-read here; Door 4 stays closed until it is.
- Door 3 needs an admission check of its own (the fixed-wage regime's criterion
  does not apply unchanged once one wage is free).
- A DELTA companion and a sectoral supply-shock design are unexplored, and
  either could be run before any kernel change.

# Revision Log

- **Version 1** (September 2026)
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- The decisive claim of Version 1 is measured: the external account closes under a tilted wage vector, at the solver-residual level (2.9e-16 at the degenerate pin, at most 7.0e-12 under a tilt), via the price-weighted clearing residual, validated against the kernel canary off equilibrium along quantity perturbations (ratio 1.0000). The reproduction table and the provenance table were updated accordingly, the open-items entry that blocked on the missing measurement was rewritten, and the code-surface note now distinguishes the closure question (answered) from the gate question (the canary must still take the vector). ADR-0021 (proposed) records the adoption decision; the reference to it was added in the reproduction section.
