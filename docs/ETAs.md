---
title: "Supply-side arm: identifying eta_s in the BETA closure"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [beyondhulten, workplan, closures, eta_s, supply-shock, identification]
last-updated: September 2026
---

This document is the workplan for the supply-side arm, the cheapest of the
candidates that restore variety inside the {BF, ALPHA, BETA} block (section
4.2 of `docs/DOCS_ASSESSMENT.md`, ranked candidate 2). House rule: plans
normally live in `ROADMAP.md` or in an ADR; this file exists at your request,
and the durable records remain the design file, the `registry/scenarios.csv`
rows, the run manifests and the ADR that authorises the `run.jl` change. The
short answer to the question that prompted it: **yes, under a supply shock
BETA separates from ALPHA and from itself** --- the row stops being a
duplicate and becomes an identification result. The pilot already exists and
its numbers are quoted below.

# Why the demand-only matrix cannot see eta_s

The reported real-GDP metric is the household-consumption Tornqvist index, and
under a demand-only shock the price block is pinned by technology alone: the
zero-profit conditions $p_i = c_i(p, w)$ are homogeneous of degree 1 in
$(p, w)$ and contain no demand term, so they determine only the ratio
$p / w = p^{*}(A)$. With the CPI numeraire, $w = 1 / \mathrm{CPI}(p^{*})$, and
at the calibrated technology ($A = 1$, hence $p^{*} = 1$) that gives $w = 1$
and $p = 1$ for **any** demand vector. The BETA supply curve
$L^{s} = \bar{L}\, [(w/P)/(w_0/P_0)]^{\eta_s}$ then returns $\bar{L}$ for
every $\eta_s$, which is why BETA reproduces ALPHA bit-for-bit (measured:
max abs difference in quantities $1.1\times10^{-16}$, in prices $0$).

Two corollaries that the session measured and that bound the design:

- A demand shock of ten times the programme size still gives $w = 1.0000000000$
  and $\max|p - 1| = 2.66\times10^{-15}$; no demand-side re-specification can
  move the real wage.
- The numeraire is not the lever either (DE-0011): re-gauging changes the
  level of $(p, w)$, not their ratio.

The only channel that moves $w/P$ is technology, i.e. `shocks.supply_shock`.

# The pilot that already exists

`experiments/probes/probe2_identification_and_shocks.jl`, section D: a
+20 \% productivity shock in sector 1 ($A_1 = 1.2$, all other sectors $A = 1$),
solved at $\theta = \epsilon = 0.5$, $\sigma = 0.9$, $\eta = 1$, warm-started
from the stored ALPHA solution. Residual below $5.4\times10^{-13}$ in every
cell.

| row | $w$ (real wage) | $L$ | real GDP F1 | real GDP F2 | real GDP F3 |
| --- | ---: | ---: | ---: | ---: | ---: |
| ALPHA ($\eta_s = 0$) | 1.0049029745 | 1.0000000000 | 1.0055526180 | 0.9881460294 | 1.0051598783 |
| BETA $\eta_s = 0.5$ | 1.0049029745 | 1.0024484897 | 1.0086846455 | 0.9912768336 | 1.0082906826 |
| BETA $\eta_s = 2$ | 1.0049029745 | 1.0098299882 | 1.0181268156 | 1.0007153159 | 1.0177291649 |
| BETA $\eta_s = 5$ | 1.0049029745 | 1.0247564458 | 1.0372202492 | 1.0198012921 | 1.0368151411 |

Three things follow, and they are the substance of the answer:

1. **BETA moves, ALPHA does not.** Employment is pinned at exactly $1$ for
   ALPHA and rises with $\eta_s$ for BETA: $+0.24$ \%, $+0.98$ \%, $+2.48$ \% at
   $\eta_s = 0.5, 2, 5$. Under F1 the headline moves from $+0.56$ \%
   (ALPHA) to $+0.87 / +1.81 / +3.72$ \% --- a spread of more than three
   percentage points across the elasticity grid, against a spread of
   $1\times10^{-16}$ in the demand-only matrix.
2. **The elasticity is recovered exactly.** $\ln L / \ln w$ returns
   $0.500000$, $2.000000$, $5.000000$; the pilot is therefore also a
   specification test of the implementation.
3. **Prices move at last** ($\max|p - 1| = 0.196$), which is the precondition
   for the other candidates (markups, capacity) to have anything to bite on.

# Design of the arm

## Shock specification

Three candidate shocks, in ascending cost. The recommendation is to run all
three, because they answer different questions.

| id | shock | what it identifies | note |
| --- | --- | --- | --- |
| `s1` | $A_1 = 1.2$, single sector (agriculture) | the pure identification test | already piloted; surgical, but economically arbitrary |
| `prog` | $A_i = 1 + \alpha \psi_i$ on the programme's own sectors ($\alpha \in \{0.05, 0.10, 0.20\}$) | the economically relevant case: a green-investment programme bids up the wages of the sectors it buys from | the impulse is 63 \% specialised construction; this shock targets exactly the margin the manuscript argues about |
| `unif` | $A_i = 1.01$ for all $i$ | the cleanest real-wage experiment; uniform, so no composition contaminates the estimate | small magnitude: a uniform 20 \% shock would move the real wage by 25 \% |

## Cells and run ids

Design `supply_etas` in `experiments/designs/supply_etas.toml`; ids follow
`<design>-<labor>-<financing>[-<variant>]` (ADR-0002). One design per shock
id keeps the manifests readable: `supply_etas_s1`, `supply_etas_prog`,
`supply_etas_unif`.

- Labour rows: ALPHA (the $\eta_s = 0$ control) and BETA at
  $\eta_s \in \{0.5, 1, 2, 5\}$, e.g. `supply_etas_s1-BETA-F2-etas2`.
- Financing: F1, F2, F3 as in the matrix, so the arm nests the published
  results.
- Shock ladder: three magnitudes per shock id, so the recovered elasticity can
  be checked for invariance.
- Count: $5 \times 3 \times 3 = 45$ cells per shock id; the pilot suggests a
  few seconds per cell, so each design is a single batch.

## What "estimating eta_s" means here

Inside the model $\eta_s$ is a parameter, not an estimate: at a solved cell it
is recovered exactly by inversion,

$$\eta_s = \frac{\ln L}{\ln (w/P)},$$

and the pilot returns the input values to six decimals. The arm therefore
delivers two publishable quantities rather than a point estimate:

1. **The identification result**: the map $\eta_s \mapsto (w/P, L, \text{real
   GDP})$ is strictly monotone and invertible under a supply shock, and
   degenerate (constant) under a demand-only shock. That is the claim the
   referee demand R2 actually needs.
2. **The sensitivity band**: the semi-elasticity of the headline to $\eta_s$,
   reported as the spread over $\eta_s \in [0.5, 5]$ per financing closure ---
   i.e. how much of the published result rests on an elasticity the data do
   not pin down.

An empirical $\eta_s$ (from an observed $\Delta \ln L$ against
$\Delta \ln (w/P)$ pair) then enters the model as a scenario, not as a
calibration; the arm is what makes that step informative.

## Preregistered signatures

Recorded with the design (ADR-0006), deviations are findings (ADR-0004):

- $\ln L / \ln (w/P)$ reproduces the input $\eta_s$ to $1\times10^{-6}$ in every
  BETA cell, and is invariant to the shock magnitude within a shock id.
- ALPHA: $L = 1$ exactly, $w > 1$; BETA: $L$ strictly increasing in $\eta_s$.
- Monotone ordering of the real-GDP headline in $\eta_s$ within each financing
  column.
- Residual gate $10^{-6}$ (mobile); budget gate $10^{-9}$.

## Gates and open items to check, not to assume

- **Canary under $A \neq 1$.** The ADR-0010 identity
  $p \cdot \text{mktc} = S - (I + X - M) + T$ has only ever been asserted at
  $p = 1$. Zero profit still holds, so the derivation should carry, but the
  terms are now valued at moving prices: assert it on the first supply cells
  as a diagnostic and record the result (an ADR if it fails).
- **Scale determinacy.** The round-gain criterion
  $\max(A_{\text{bill}}/\lambda + (1-m)(1-s)\,fs) < 1$ uses the baseline
  intermediate coefficients; a productivity shock changes the optimal
  input mix at $\epsilon = 0.5$, so the criterion must be re-evaluated for the
  arm.
- **BF in the arm.** The recombination corner ($\eta = 0$ with
  $\eta_s > 0$) becomes informative only here, but it inherits the
  identification defect of ADR-0017: report its aggregates, never its
  quantities, until sector-specific wages land.

# Workplan

1. **ADR-0018 (blocking).** `experiments/run.jl` hard-codes
   `Shocks(ones(N), ones(N), zeros(N))` in `build_cell_model` and
   `solve_cell`; the arm needs a per-cell supply-shock source and magnitude
   from the design file. The ADR fixes the schema, the scenario-id convention
   for $\eta_s$ variants, and the rule that the baseline reference stays the
   no-shock one (real GDP remains "against the baseline").
2. **Bridge, before the ADR lands.** Extend the pilot to the full grid as
   `experiments/probes/probe5_supply_arm.jl` (read-only, no `src/` change), so
   the numbers exist and the signatures can be preregistered against
   measurement rather than expectation.
3. **Registry.** Append `planned` rows to `registry/scenarios.csv` for every
   cell (45 per shock id), with the parameter cells, the shock id and
   magnitude, and `eta_s`.
4. **Design and preregistration.** Write
   `experiments/designs/supply_etas_<id>.toml`, then

   ```{.bash}
   julia --project=. experiments/run.jl --list supply_etas_s1
   julia --project=. experiments/run.jl --preregister supply_etas_s1 --actor hermes-agent
   ```

   and commit the preregistration record before any cell is created.
5. **Execution.**

   ```{.bash}
   julia --project=. experiments/run.jl --design supply_etas_s1 --actor hermes-agent
   ```

   then `supply_etas_prog` and `supply_etas_unif`.
6. **Gates after every batch.**

   ```{.bash}
   julia --project=. scripts/status.jl        # 0 warnings
   julia --project=. scripts/check_repo.jl    # 0 violations
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```
7. **Write-up.** A paper table citing the `supply_etas-*` run ids: rows
   ALPHA and BETA ($\eta_s$ grid), columns F1/F2/F3, cells = real GDP and
   employment, plus the recovered elasticity column. Update
   `docs/DOCS_ASSESSMENT.md` (section 4.2, the ranked proposal) and
   `registry/closures.toml [labor.BETA]` open gates; one
   `docs/log/2026-09.md` entry.

# What the arm does not do

The arm moves prices **through technology only**: the zero-profit block
contains no demand term, so $p \neq 1$ requires $A \neq 1$. Put the
demand-only programme back in and the matrix returns to $p = 1$, $w = 1$,
$L = \bar{L}$ --- BETA collapses onto ALPHA again. The supply-side arm is
therefore a **different experiment**, not a repair of the published 5 x 3
matrix: it identifies $\eta_s$ and it prices the elasticity, but it does not
restore variety inside the demand-shock block.

Variety *in the matrix the paper publishes* requires the price block to depend
on demand, i.e. one of:

- **demand-sensitive factor prices** --- sector-specific wages (ranked
  candidate 1): N sectoral labour markets make $w_i$, and through zero profit
  $p_i$, functions of sectoral labour demand. Zero profit is preserved (no
  rents to dispose of), it repairs the ADR-0017 defect, and it is the only
  variant that is *earned* by the model rather than assumed: the new
  parameters are sectoral supply elasticities of the same family as $\eta_s$.
- **a demand-sensitive markup** $\mu_i = 1 + \kappa\,(y_i/y_{i0} - 1)$:
  nested at $\kappa = 0$ (the current model), a few lines of code, and a clean
  sensitivity ladder. But markup income is profit, so it needs a
  profit-income-and-spending closure (who receives the rent, and how it is
  spent) or the accounting identity behind the canary breaks --- that is the
  real cost, not the code.
- **a fixed factor with capacity**: the standard short-run mechanism, but the
  largest surgery --- a second factor, a rental schedule and a rental-income
  closure.

Once prices respond to demand, BETA separates from ALPHA in the demand matrix
as well: the numeraire then delivers a real wage that moves with the shock, so
$L^{s} = \bar{L}\,(w/P)^{\eta_s}$ bites. The sequencing this implies:
supply-side arm first (cheap, earns the identification claim for $\eta_s$),
then sector-specific wages (the repair that gives the published matrix
variety), with the markup ladder kept as a robustness arm rather than as the
headline mechanism.

# Risks

- **Solver.** The pilot used `solve_beta` with eight continuation rungs and
  passed at $\eta_s = 5$; the `prog` shock at $\alpha = 0.20$ hits the
  construction block harder, so budget more rungs and keep the warm start.
- **Magnitude.** A uniform shock must stay small (one to five per cent); the
  single-sector shock can be large because its consumption weight is small.
- **Interpretation.** The arm changes the *experiment*, not the closure;
  nothing in `src/` moves. It does not fix the BF row (that is candidate 1,
  sector-specific wages), and it does not make $\eta_s$ identified without a
  supply shock --- the demand-only matrix stays a null result for this
  parameter, which is itself the finding.

# Deliverables

- `experiments/designs/supply_etas_{s1,prog,unif}.toml` with preregistration
  records.
- One run manifest per cell under `runs/supply_etas-*/`, rows in
  `registry/scenarios.csv`, index rows in `runs/index.csv`.
- ADR-0018 (design schema for supply shocks) and, if the canary or the
  round-gain check fails under $A \neq 1$, a follow-up ADR.
- `paper/tables/supply_etas_flows.md`, generated from the manifests and citing
  the run ids (ADR-0004).
