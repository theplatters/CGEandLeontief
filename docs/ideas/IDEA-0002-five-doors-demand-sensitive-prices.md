---
title: "IDEA-0002: the five doors to demand-sensitive prices (with the Verdoorn arm)"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, prices, demand-sensitivity, verdoorn, five-doors, identifiability]
last-updated: September 2026
---

**Version 1** (September 2026)

Status: `unexplored`. The menu is derived and ranked in `docs/ETAs.md`; nothing
in it has been probed, calibrated or run.

# The shared derivation

The doors are one family, not five unrelated ideas, because they all answer the
same obstruction. The zero-profit block is homogeneous of degree one in $(p, w)$
and contains no demand term, so under $A = 1$ every regime sits at
$p = w = \mathrm{CPI} = 1$ and the real wage is pinned at its anchor; the
quantity block then reduces to an affine multiplier system in which none of
$(\theta, \epsilon, \sigma, \eta_s)$ appears. That is Lemma 1 and Corollary 1 of
`paper/equivalence.tex` (the price block is demand-invariant; the real wage is
invariant), and it is why ALPHA is equivalent to BETA and GAMMA to DELTA as
propositions rather than as numerical coincidences.

Variety inside the demand-only matrix therefore requires **demand to enter the
price block**. Every door below does that in a different place, and each carries
its own price sign, its own parameter and its own accounting cost.

# The doors

- **Labour: sectoral wages.** $N$ sectoral labour markets make $w_i$, and through
  zero profit $p_i$, functions of sectoral labour demand. Zero profit is
  preserved, no rents are created, and the new parameters are sectoral supply
  elasticities of the same family as $\eta_s$ --- the only door whose parameters
  are *earned by the model rather than assumed*. Its all-rigid corner
  ($\eta_{s,i} = 0$ for every sector, so $L^{cm}_i = \bar L_i$) is exactly the
  `eta = 0` row of the current matrix, which makes the existing BF cells a free
  endpoint for the family. `docs/ETAs.md` ranks it candidate 1 for the published
  matrix, and `docs/WORKPLAN_SENSITIVE_PRICES.md` is the implementation route.
- **Markup: a demand-sensitive margin.** $\mu_i = 1 + \kappa (y_i / y_{i0} - 1)$,
  nested at $\kappa = 0$. A few lines of code, but markup income is profit, so it
  needs a profit-income-and-spending closure or the accounting identity behind
  the canary breaks. The cost is the closure, not the code.
- **Capacity: a fixed factor.** The standard short-run mechanism and the largest
  surgery: a second factor, a rental schedule, a rental-income closure.
- **Capacity's reduced form: an endogenous utilization externality.** The
  effective productivity $A^{\mathrm{eff}}_i = A_i (y_i / \lambda_i)^{-\delta}$
  inside the unit-cost hook. Demand raises marginal cost, zero profit is
  preserved and **no rent is created**, so there is no income-closure problem at
  all: one line in the cost function, one new parameter. `docs/ETAs.md` calls it
  the cheapest door of all.
- **External: endogenous import prices.** $p^M_i = p^M_{i0} (M_i / M_{i0})^{\zeta}$
  with $\zeta \geq 0$, an upward-sloping foreign supply curve. Import prices enter
  the intermediate index and the CPI through the margins, so relative prices
  become functions of import volumes and hence of demand. Moderate surgery, and
  it exploits precisely the external-account machinery of ADR-0019: the
  programme's import content bids up import prices, a structuralist absorption
  and terms-of-trade channel.

Not a door: supply-side technology shocks ($A \neq 1$). They separate BETA from
ALPHA but do not repair the demand-only matrix, and they are recorded separately
as `docs/ideas/IDEA-0001-climate-productivity-shocks.md`. Also not a door, by
proof: any rule that writes the *nominal* wage as a function of endogenous
aggregates, including nominal indexation, cannot create demand-sensitive real
outcomes (Corollary 3 of `paper/equivalence.tex`, the no-go for nominal wage
rules). That family is closed.

# The Verdoorn arm in particular

The utilization externality is the door worth looking at first, for four reasons.

- **It is the cheapest.** One term in the unit-cost hook, one parameter, no new
  equation, no new closure, no rent to dispose of.
- **It has a sign.** With $\delta > 0$ demand raises marginal cost, the
  utilization or congestion reading. With $\delta < 0$ it is the
  Kaldor--Verdoorn case: increasing returns make prices *fall* with demand. No
  other door gives a negative price response, which makes the arm a genuine
  bracket rather than a robustness check: the programme's price effect becomes a
  two-channel question, pass-through against productivity, instead of a
  one-signed prediction.
- **It is honest about its status.** `docs/ETAs.md` calls it a reduced form. The
  parameter is *assumed*, not derived from a sectoral supply elasticity, so it
  cannot carry the identification claim that the labour door carries; it is a
  sensitivity ladder, and the ladder is the point.
- **It composes.** The utilization term and the sectoral-wage term can be active
  together, so the arm can be run as a companion to the labour door rather than
  as a substitute, and the two price responses can be separated by switching one
  off.

What it would need: the same per-cell parameter source in the design schema that
`docs/ETAs.md` records as the blocking item of its own workplan (the harness
hard-codes a null shock and exposes no per-cell parameter), a calibration
argument for $\delta$ that this note does not supply, and the same diagnostic as
every other door: does `max abs(p - 1)` move with the financing cell.

# Why one note and not five

The doors share the derivation, the diagnostic and the selection criterion, and
they are complements rather than rivals: the labour door supplies the earned
parameters, the Verdoorn arm supplies the opposite sign, the external door
supplies the open-economy channel that the ADR-0019 account was built for, and
the markup and fixed-factor doors supply robustness. Splitting them into five
notes would hide the fact that choosing among them is a single decision about
where demand is allowed to enter the price block.

The criterion for choosing, when the time comes, is the one this repository
already uses: what the model earns, what the accounting absorbs, and what the
paper can defend to a referee. By that criterion the labour door leads, the
Verdoorn arm is its cheap companion, the external door is the natural
open-economy extension, and the markup and fixed-factor doors wait.

# Revision Log

- **Version 1** (September 2026)