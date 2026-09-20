---
title: "Supply-side arm: identifying eta_s in the BETA closure"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [beyondhulten, workplan, closures, eta_s, supply-shock, identification]
last-updated: September 2026
---

**Version 3** \textcolor{revisionV2}{(September 2026)}
**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

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

# Status of the arm \textcolor{revisionV2}{\normalsize [v3, 2026-09-20]}

\textcolor{revisionV2}{Version 3 records three things: the schema ADR is
drafted (ADR-0023, proposed); the bridge step is done --- the pilot is
re-measured on the sectoral kernel by \texttt{probe15} (2026-09-19), whose
numbers supersede the table above as the current evidence; and an
identification-vehicle decision is fixed. The arm still changes the
experiment, not \texttt{src/}: no closure, no kernel, no generation.}

\textcolor{revisionV2}{\textbf{The bridge measurements (probe15, full-71,
$A_1 = 1.2$; recorded in \texttt{docs/log/2026-09.md}).} The $\eta_s$ rungs
separate: employment 1.000000 / 1.002187 / 1.008218 at $\eta_s$ = 0 / 0.5 / 2
(sectoral uniform vector, F2), against a 0.02\,\% employment movement over
the same ladder under the demand-only programme; the real-wage spread reads
1.215 / 1.105 / 1.043. The rigid-corner nesting (Proposition 3) is
\emph{shock-independent}: the canary is bit-exact with and without the
shock. Two cautionary measurements: the shock alone moves prices
($\max|p-1| = 0.226$ with no programme), so under a supply shock
$\max|p-1|$ is not a demand-sensitivity measure --- the demand signature
remains the F1 vs F2 = F3 difference; and the welfare effect flips sign
across the rows (supply-only +1.20\,\%, supply + programme -0.40\,\%,
demand-only -1.57\,\%), i.e. the two channels do not add linearly.}

\textcolor{revisionV2}{\textbf{Identification vehicles (the v3 decision).}
The primary identification claim runs on the \textbf{scalar single-wage
BETA} cells (2N+2, the executed matrix form), because the preregistered
signature is scalar: $\ln L / \ln (w/P)$ recovers the input $\eta_s$ to
$10^{-6}$ in every cell and is invariant to the shock magnitude within a
shock id. Each design additionally carries the \textbf{sectoral uniform
vector} cells (3N+1, $\eta_{s,i} = \eta_s$ for every $i$) as the robustness
row: there the identification is a vector ($\eta_{s,i} =
\ln(L_i/\bar L_i)/\ln(w_i/\Pi)$ per sector) and no scalar signature applies;
report per-sector recovery, never one number. Both forms are live in the
kernel and in the executed v9 design (matrix cells = scalar, ladder cells =
sectoral), so no promotion is needed. The ADR-0023 schema records this and
keeps the two rows separate in every table.}

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
  \textcolor{revisionV2}{Status (v3): partly measured, not assumed ---
  probe15 row set 3 asserts the rigid-corner nesting (Proposition 3)
  bit-exact with and without $A_1 = 1.2$, and the external-account identity
  is asserted on the first arm cells as planned.}
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
   \textcolor{revisionV2}{Status (v3): the schema ADR is drafted as
   \texttt{ADR-0023} (proposed, 2026-09-20) --- per-design \texttt{[shock]}
   block, default \texttt{kind = "none"} (existing designs bit-identical),
   the ADR-0002 id convention with the shock id in the variant slot, and the
   no-shock reference unchanged.}
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

# What the arm does not do \textcolor{revisionV1}{\normalsize [extended v2]}

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
- \textcolor{revisionV1}{\textbf{the capacity door's reduced form} --- an
  endogenous utilization externality $A^{\mathrm{eff}}_i = A_i\,(y_i /
  \lambda_i)^{-\delta}$ inside the unit-cost hook: demand raises marginal
  cost, zero profit is preserved, and \emph{no rent is created}, so there is
  no income-closure problem at all --- a one-line change in the cost function
  with one new parameter. With the opposite sign ($\delta < 0$) it is the
  Kaldor--Verdoorn case: increasing returns make prices \emph{fall} with
  demand --- variety in the opposite direction. The cheapest door of all.}
- \textcolor{revisionV1}{\textbf{the external door} --- endogenous import
  prices $p^M_i = p^M_{i0}\,(M_i / M_{i0})^{\zeta}$ with $\zeta \geq 0$: an
  upward-sloping foreign supply curve. Import prices enter the intermediate
  index and the CPI through the margins, so relative prices become functions
  of import volumes and hence of demand. Moderate surgery, and it exploits
  precisely the external-account machinery of ADR-0019: the programme's
  import content bids up import prices --- a structuralist absorption /
  terms-of-trade channel. New candidate.}
- \textcolor{revisionV1}{\textbf{the technology door} --- supply shocks
  $A \neq 1$: this arm, the present workplan; it separates BETA from ALPHA
  but does not repair the demand-only matrix.}

\textcolor{revisionV1}{\textbf{The five doors (added v2).} The ex ante
derivation now lives in \texttt{paper/equivalence.tex}: the zero-profit
block is homogeneous of degree one in $(p, w)$ and demand-free (Lemma 1
there), so under $A = 1$ every regime sits at $p = w = \mathrm{CPI} = 1$
and the real wage is pinned at its anchor --- the quantity block reduces to
one affine multiplier system in which none of $(\theta, \epsilon, \sigma,
\eta_s)$ appears. Variety therefore requires demand to enter the price
block, and the list above is exactly the door menu: labour, markup,
capacity (with its reduced form), external, technology.}

\textcolor{revisionV1}{\textbf{Ruled out ex ante} (the no-go corollary,
Corollary 3 of \texttt{paper/equivalence.tex}): no rule that writes the
\emph{nominal} wage as a function of endogenous aggregates --- an aggregate
wage curve $w = \bar{w}\,(\alpha^{\!\top} y / \bar{L})^{\varphi}$, nominal
indexation, or a numeraire change --- can create demand-sensitive real
outcomes. The price block still forces $p = w\,\pi(A)$ and
$\mathrm{CPI}(p) = w\,\mathrm{CPI}(\pi)$, so the real wage is invariant to
whatever sets $w$; in the fixed-wage gauge such a rule is either vacuous or
over-determines the system (it forces $\alpha^{\!\top} y = \bar{L}$,
silently reverting GAMMA to ALPHA). A rule on the \emph{real} wage is the
BETA supply equation rearranged and is already in the model, where the
demand-only design pins it to the anchor. Quantity-side levers (endogenous
saving rates, endogenous export demand) move the multiplier, not the price
block, and add no closure variety. Do not schedule this family.}

Once prices respond to demand, BETA separates from ALPHA in the demand matrix
as well: the numeraire then delivers a real wage that moves with the shock, so
$L^{s} = \bar{L}\,(w/P)^{\eta_s}$ bites. The sequencing this implies:
supply-side arm first (cheap, earns the identification claim for $\eta_s$),
then sector-specific wages (the repair that gives the published matrix
variety), with the markup ladder kept as a robustness arm rather than as the
headline mechanism.

\textcolor{revisionV1}{\textbf{Sequencing update (v2).} The ranking now
reads: the supply-side arm first (unchanged --- cheap, earns the
identification claim for $\eta_s$); then sector-specific wages (the labour
door, still candidate 1 for the published matrix --- it repairs the BF
row's common-wage defect); the utilization/Verdoorn externality as the
cheap third arm (no income-closure problem; with $\delta < 0$ it delivers
Kaldorian variety in the opposite direction); the terms-of-trade arm as the
open-economy candidate with the best fit to the ADR-0019 external account;
the markup ladder and the explicit fixed factor kept as robustness arms.
The nominal-wage-rule family is ruled out ex ante (see above) and stays off
the workplan.}

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

# Revision Log

- **Version 1** (September 2026)
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Section 5
  extended: the three demand-sensitive-price candidates completed to the
  five-door menu (new: the external door --- endogenous import prices /
  terms of trade; the capacity door's reduced form --- the
  utilization/Verdoorn externality with no income-closure problem); the
  nominal-wage-rule family ruled out ex ante (Corollary 3 of
  \texttt{paper/equivalence.tex}); sequencing updated.
- **Version 3** \textcolor{revisionV2}{(September 2026)} --- The schema ADR
  is drafted as ADR-0023 (proposed); the bridge step is done: probe15
  re-measures the pilot on the sectoral kernel (rung separation, the
  shock-independent nesting, the price-response and welfare-sign
  findings); the identification-vehicle decision is fixed (scalar
  single-wage BETA carries the scalar signature, the sectoral uniform
  vector the per-sector robustness row); workplan step 1 and the canary
  gate updated.
