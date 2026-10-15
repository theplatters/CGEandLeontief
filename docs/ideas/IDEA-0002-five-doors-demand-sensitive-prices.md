---
title: "IDEA-0002: the five doors to demand-sensitive prices (with the Verdoorn arm)"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, prices, demand-sensitivity, verdoorn, five-doors, identifiability]
last-updated: October 2026
---

**Version 3** (October 2026)
**Version 2** (October 2026)
**Version 1** (September 2026)

Status: `promoted`. The capacity and markup doors are adopted by ADR-0029 and
ADR-0030 as the pricing closures `P1` and `P2` (`registry/closures.toml`) and
executed in the `pricing_*_v1` designs; cite those run ids, not this note. The
probe numbers in "Measured" below remain uncitable. The external (import-price)
door, the fixed-factor door and the technology door are not promoted and stay
ideas inside this note.

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
  matrix, and `docs/VariationinGamma.md` (the wage-structure door) plus `docs/ETAs.md` are the implementation routes.
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

# Measured: the capacity and markup doors (probe20, October 2026)

`experiments/probes/probe20_pricing_door.jl` measures ROADMAP Phase 7 before
any ADR, on the full-71 A-bill calibration, the 2024 programme
($G_0 = 0.01331$), $\theta = \epsilon = 0.5$, $\sigma = 0.9$, all three
financing columns. These are probe numbers, not run ids.

- **Route 1, capacity / Kaldor--Verdoorn.** $A^{\mathrm{eff}} = (y/\lambda)^{-\delta}$ is solved
  as probe17's damped fixed point, with $\delta$ from $-0.75$ to $0.75$ in
  steps of $0.125$.
- **Route 2 Stage A, retained markup.** $p = \mu c$ with
  $\mu = 1 + \kappa (y/\lambda - 1)$ and the household budget unchanged, for
  $\kappa$ from $0$ to $0.75$. The markup needs **no kernel change** to
  measure. Every first-order-condition term in the kernel's demand hook
  carries $p_u^{\epsilon} y_u$, so evaluating the hook at $y/\mu^{\epsilon}$
  puts marginal cost $c = p/\mu$ into the first-order conditions, while
  purchaser prices stay in the price indices. This is an exact form of the
  ROADMAP's "FOC correction". It must be applied at every site that carries
  $p^{\epsilon}$: the demand hook, the external-account canary and the GDP
  components.

**Canaries.**

- At $\delta = \kappa = 0$, the fixed gauge reproduces
  `matrix_5x3-v10-GAMMA-*` (employment 0.999027 / 1.001260 / 1.017796) and
  the mobile gauge reproduces `matrix_5x3-v10-ALPHA-*` (consumption
  +0.143 % / -1.823 % / -1.823 %).
- Route 1 is rent-free: the external-account gap stays below $4 \cdot 10^{-11}$.
- Route 2's gap equals $-\Pi_{ret}$, with
  $\Pi_{ret} = \sum_i (\mu_i - 1) c_i y_i$, to $1.1 \cdot 10^{-11}$ on every
  cell in both gauges. The retained profit is therefore the only missing term
  in the identity, and this is the sign it is booked with.
- All 75 cells solve.

**Findings.**

1. **The sign of the price response depends on how the programme is
   financed.** Kaldor--Verdoorn ($\delta < 0$) lowers the GDP deflator only
   under F3: 0.9853 at $\delta = -0.5$, against 1.0035 (F1) and 1.0013 (F2).
   Under budget-neutral financing, the programme shrinks the sectors it
   displaces. Under increasing returns those sectors lose productivity, and
   their price rise outweighs the fall in the programme sectors. At
   $\delta = -0.75$ the deflator is 1.0099 (F1) and 1.0062 (F2). The bracket
   ("raises or lowers prices depending on the door") is therefore conditional
   on the financing closure, which strengthens "which money?" and qualifies
   the planned 6.4 sentence.
2. **$\max|p - 1|$ cannot carry the bracket.** It is nearly symmetric in
   $\delta$ (F3: 0.1280 at $-0.5$, 0.1285 at $+0.5$), and small sectors
   (glass, ceramics) dominate it. The Phase 7 table needs a signed aggregate,
   either the deflator or the $\lambda$-weighted mean price change.
3. **Same prices, different quantities.** Under a Hicks-neutral $A$, the two
   doors scale unit cost by $(y/\lambda)^{\delta}$ and
   $1 + \kappa(y/\lambda - 1)$. These coincide to first order at
   $\kappa = \delta$, and the sectoral price patterns are nearly identical.
   The doors differ in where the wedge goes. Capacity absorbs it as extra
   input use; markup absorbs it as a rent that leaks from the household
   circuit. At magnitude 0.5 under F3:

   | | Deflator | Employment | Consumption | $\Pi_{ret}/Y$ |
   | --- | --- | --- | --- | --- |
   | GAMMA ($\delta = \kappa = 0$) | 1.0000 | 1.0178 | +2.26 % | 0 |
   | Capacity $\delta = +0.5$ | 1.0215 | 1.0388 | +2.20 % | 0 |
   | Markup $\kappa = +0.5$ | 1.0095 | 1.0108 | +0.41 % | 0.95 % |

   The rent/no-rent contrast planned for 6.4 is therefore a contrast in
   quantities and welfare, not in prices. Under retention, the markup's
   welfare cost is largely hoarding, which is the case for Stage B.
4. **Cobb--Douglas pass-through canary.** At $\theta = \epsilon = 1$ the price
   block is exactly $\log p = (I - B)^{-1} r$, with
   $B = \mathrm{diag}(1 - fs)\,\Omega$ and $r = \log \mu$ (markup) or
   $\delta \log(y/\lambda)$ (capacity). It holds to $2.5 \cdot 10^{-12}$ on
   every probed cell. A constant uniform markup gives the closed form
   $p = \bar\mu^{(I-B)^{-1}\mathbf 1}$ to $7 \cdot 10^{-16}$. The cascade
   $(I-B)^{-1}\mathbf 1$ ranges from 1.42 to 4.36 (its $\lambda$-mean is
   2.18), so a uniform 5 % markup raises prices by 11.3 % on the
   $\lambda$-weighted mean.

   **Correction to the ROADMAP's Route 2 canary 6:** "a uniform $\kappa$
   scales prices by exactly $1 + \kappa(y/y_0 - 1)$" ignores the
   intermediate cascade. It misses by 0.164 in $\log p$ at
   $\bar\mu = 1.05$. The canary should be stated as the Cobb--Douglas closed
   form above.
5. **Numerics.** On the CES ladder, iterations rise from about 30 near
   $\delta = 0$ to 70--76 at $\delta = -0.5$ and 158--181 at $\delta = -0.75$.
   No rung failed. Under Cobb--Douglas the $\delta = -0.5$ iteration
   *diverges* at damp 0.6 and converges at damp 0.2. The equilibrium exists,
   so the failure belongs to the iteration, not the model. A failed fixed
   point at one damping is therefore not evidence of non-existence, and the
   uniqueness guard (multi-start plus the ADR-0017 check at the converged
   $A$) remains unmeasured.
6. **Mobile gauge (ALPHA, $\eta = 1$, CPI = 1).**
   - F2 = F3 neutrality survives both doors (identical rows).
   - Prices move relatively: $\max|p-1|$ is 0.10--0.13.
   - The deflator falls at $\delta = -0.5$ in every column (0.997) and rises
     at $\delta = +0.5$ and $\kappa = 0.5$ (1.003).
   - The real wage rises slightly under every door (+0.08 % to +0.20 %).
     That sign is measured but not yet interpreted.
   - The markup's consumption cost is small here (F2: -1.90 % against
     -1.82 % at $\kappa = 0$), because employment is pinned.

**Not yet measured:** the $\delta < 0$ uniqueness guard; Stage B (domestic
profit, the $s_\pi$ nesting canary); and an independent check of
$w_i = c_i\,\partial\Phi_i/\partial L_i$. The emulation satisfies that last
identity by construction, so a real check needs the kernel hook.

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
- **Version 3** (October 2026) --- status `scoped` to `promoted` for the capacity and markup doors (ADR-0029, ADR-0030; closures P1/P2; designs `pricing_{capacity,markup,profit,mobile}_v1`); the external, fixed-factor and technology doors remain ideas.
- **Version 2** (October 2026) --- status `unexplored` to `scoped` (ROADMAP Phase 7 specifies the markup and capacity doors); added "Measured" (probe20: both Phase 7 routes on F1/F2/F3, the Cobb--Douglas pass-through canary and the mobile-gauge row), including the financing-dependent sign of the Kaldor--Verdoorn price response and the correction to the ROADMAP's uniform-$\kappa$ canary.
