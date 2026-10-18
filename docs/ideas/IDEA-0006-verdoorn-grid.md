---
title: "IDEA-0006: the Verdoorn grid --- which matrix closures to carry the capacity door"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-10-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, pricing, verdoorn, dissolution, equivalence, closures, phase-7]
last-updated: October 2026
---

**Version 1** (October 2026)

Status: `scoped`.

# The question

Which of the five labour-market closures of the 5x3 matrix should carry the
capacity door (closure `P1`, the Kaldor-Verdoorn / utilisation elasticity), given
that it is at present executed in two gauges only. The question is not how to
compute more cells; it is which extra cells retire a *claim* and which only add
numbers.

# What exists

Measured from `runs/` (2026-10-18). Sixty-nine pricing-door runs, all of them in
two gauges:

| Labour closure | Runs | Composition |
| --- | --- | --- |
| GAMMA (fixed real wage) | 45 | capacity v1 9, capacity v2 15, markup 9, profit 12 |
| ALPHA (mobile labour) | 24 | mobile v1 12, mobile v2 12 |
| BETA, BF, DELTA | 0 | --- |

Every run carries the uniqueness gate (multi-start spread, full-rank residual
Jacobian); the nine nesting cells (`delta = 0`, `kappa = 0`) are exactly those
that reproduce the demand-only matrix, and sixty runs have moving prices.
`registry/closures.toml` `[pricing.P1]` already states that the door "dissolves
the GAMMA = DELTA and ALPHA = BETA equivalences of Lemma 1 outside the
price-invariant-technology slice" --- claimed in the registry, asserted in
Section 6.4's closing paragraph, and measured in neither pair. The no-reallocation
equivalence is the third (Proposition `nest`).

# The three equivalences and their door counterparts

| Equivalence (demand-only only) | What differs | Counterpart cell | Verdict |
| --- | --- | --- | --- |
| GAMMA = DELTA (Prop. `gd`) | substitution elasticities at the Leontief limit | DELTA x door | thin, see below |
| ALPHA = BETA (Prop. `ab`) | one aggregate real-wage supply equation vs N sectoral ones | BETA x door | the one to run |
| BF = sectoral ladder endpoint (Prop. `nest`) | frozen allocation | BF x door | optional, expensive |

# Why GAMMA and ALPHA came first and stay the default

They are the two poles of the margin the door acts on. In GAMMA the real wage is
pinned, so a demand-driven change in unit cost must appear in prices --- the
door's cleanest display, and the gauge whose demand-only counterpart reports no
price movement at all. In ALPHA prices are pinned by the CPI numeraire, so the
same change must appear in quantities, which is the gauge that answers whether
financing neutrality survives the door. They are also the two gauges that nest
the demand-only matrix exactly (`delta = 0` reproduces
`matrix_5x3-v10-GAMMA-*` and the mobile rows bit for bit), which makes them the
cheap canaries.

# Why DELTA is thin under the door

GAMMA and DELTA share the same labour condition, `w / Pi(p) = wbar = 1`; DELTA adds
a *production-side* restriction, the Leontief limit, executed as
`theta = epsilon = sigma = 1e-4` with `delta_epsilon = 1e-4` (see the
`matrix_5x3-v10-DELTA-*` cells). On the demand-only slice that restriction is
invisible --- the quantity block is elasticity-free (Lemma 2), which is exactly
why the coincidence with GAMMA is exact rather than approximate. The door is what
switches the restriction back on, so a DELTA cell under the door is not a
labour-closure experiment at all: it measures the substitution axis, at one
point, with the wage pinned. That axis is one this revision deliberately does not
vary anywhere (the elasticity-gradient material went with the old draft), so the
cell invites the fair question of where the rest of the gradient is. Its value is
therefore narrow: it closes Proposition `gd`'s dissolution sentence and nothing
more.

# Why BETA is the one to add

BETA differs from ALPHA by a genuine *matrix* distinction --- one aggregate
real-wage supply equation instead of N sectoral ones (the BETA cells carry
`eta_s = 0.5`, 2N+2 unknowns) --- which is the axis Section 5 built the matrix on.
Section 6.4 asserts this dissolution in prose ("elastic labour supply differs from
the standard CGE closure only through the real wage, which the demand-only matrix
pins; a pricing door releases it"). The cheapest measured counterexample is one
`delta` value across the three financing columns, i.e. three cells.

# Why not BF, and how much of that is computational

The conceptual half is the stronger one. BF is the column where prices *already*
move without any door (the deflator by +0.63 / +0.72 percent, sectoral dispersion
up to 0.22 percent). The door's whole purpose in Section 6.4 is to make prices
move where they otherwise cannot, so a BF cell measures the interaction of two
open price channels rather than new variety. Its equivalence partner is also not a
headline closure but the sectoral ladder endpoint (Prop. `nest`), and the ladder
cells exist on the demand-only slice only.

The computational half is real: at `eta = 0` the system is the largest in the
kernel (unknowns `[p; y; w(1:N); F]`), the registry calls BF-F1 "the stiffest
cell", and the `eta = 0` external booking is the one that needed ADR-0020 option
C's redefinition. With `delta < 0` a positive feedback, BF x door is precisely
the combination in which the uniqueness gate is likeliest to fail --- and the
registry treats non-convergence at extreme `delta` as a documented finding, so a
failure there would have to be written up rather than reported.

# Proposed specification

A new design, `experiments/designs/pricing_dissolution_v1.toml`, on the
`pricing_capacity_v2.toml` schema (`labor` is a per-cell field, so no kernel work
is expected --- to be confirmed on the smoke run).

- Cells: BETA x {`delta = -0.5` goods, `delta = +0.5` goods} x {F1, F2, F3} = 6; DELTA x the same 6, with the Leontief limits copied from `matrix_5x3-v10-DELTA-*`; optional BF x {`delta = -0.5` goods} x {F1, F2, F3} = 3.
- Nesting canaries: `delta = 0` must reproduce `matrix_5x3-v10-{BETA,DELTA,BF}-F*` bit for bit.
- Gates: the matrix gates plus the uniqueness gate on every priced cell; continuation in `delta` (`continuation_steps = 4`).
- Cost: 12 cells (15 with BF), about a sixth of the existing 69, plus one ADR and one registry entry for the new generation.
- Untested combination to watch: DELTA at the 1e-4 limit with `delta < 0`. A failure there is a finding about the limit, not about the closure.

# What it buys, and what it does not

It converts two asserted dissolutions into one measured counterexample each, in an
executed generation, so that Section 6.4's scope sentence ("the equivalences hold
on the slice `A = 1` *and* price-invariant technology") cites cells instead of an
argument. It does not extend the variety claim to the full 5x3 grid and must not
be written as if it did; no new closure is introduced, the door stays rent-free,
and no income-closure work follows from it.

# Adjacent item found while counting

Appendix J's "All 42 priced cells pass" is not reproducible from the manifests:
there are 69 door runs, 60 with moving prices, 45 in the GAMMA gauge, 9 nesting.
42 is consistent with the GAMMA gauge minus nesting (45 minus 3), but a reader
cannot reconstruct that. Either name the set or use the countable statement.

# Open decisions for the operator

- File and run the 12-cell grid, or only the 6 BETA cells, or neither (in which case Section 6.4 and Appendix E get the one-sentence scope declaration instead).
- Include BF x 3, or leave BF out with its reason stated in the caption.
- Whether the dissolution cells are reported in Appendix E or only cited by run id.
