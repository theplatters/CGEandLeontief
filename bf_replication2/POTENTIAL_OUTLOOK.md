---
title: "Potential Outlook: Strategies for the Remaining HtM Convergence (NaN) Problems"
project: "beyond-hulten / bf_replication2"
date: "August 2026"
tags: [replication, solver, homotopy, nlsolve, path, htm, outlook]
---

# Current State of the NaN Problem

The HtM sweep (loop 2) is fully populated, but two of the six HtM shares did not
reach the terminal shock scale $t = 1.0$ and were **reconstructed** from the last
convergent $t$-point instead of genuinely converged:

| $\phi_{HtM}$ | regime | last convergent $t$ | status |
|---|---|---|---|
| 0.2 | Cobb-Douglas | 0.70 | reconstructed (underestimate) |
| 0.8 | benchmark | 0.99 | reconstructed |

The reconstruction is a reporting stopgap, not a solution. The CD cell at
$\phi_{HtM} = 0.2$ is known to be an underestimate of the true $t = 1.0$
magnitude. The goal below is to convert these two cells into real converged
solutions so the "Reconstructed" flags in `VALIDATION.md` can be removed.

See `VALIDATION.md` (HtM sweep section) and `archive/STATUS.md` for the current
numbers and the reconstruction disclaimer.

# Root Cause: Solver Deviation (Route A vs Route B)

`archive/WORKPLAN.md` specified **JuMP + PATHSolver (Route B)** as the intended
solver. The implemented `src/model.jl` instead uses **NLsolve +
Fischer-Burmeister (Route A fallback)** because the container lacked the memory
for the full finite-difference Jacobian. This deviation is acknowledged in
`archive/STATUS.md` and `VALIDATION.md`.

Why it matters for the NaN cells: PATHSolver is a dedicated mixed-complementarity
solver and is substantially more robust near the sticky-labor complementarity
($keynes = -1$). The hand-rolled NLsolve + FB reformulation is the most likely
source of the two non-convergences. The solver choice should therefore be the
first suspect when designing a fix.

# Candidate Strategies

Each option below lists the idea, its pros and cons, a rough effort estimate,
and how directly it targets the two failing cells.

## Route B Restoration: JuMP + PATHSolver

- Idea: replace NLsolve + FB with the originally planned MCP solver by wrapping
  the equilibrium system in JuMP's `MathOptInterface` / MCP form and calling
  PATHSolver.
- Pros: principled; purpose-built for complementarities; most likely to resolve
  both failing cells directly.
- Cons: requires greater than 15 GB RAM for the 668-variable FD Jacobian (run
  on the host Mac, not the container); needs a JuMP/MCP reformulation of the
  model.
- Effort: medium-high.
- Applicability: highest-confidence fix; the natural baseline target.

## Homotopy Continuation (enhanced)

- Idea: deform the problem from a solved simple case to the target, using $t$
  as the homotopy parameter. A basic version already exists as the continuation
  refinement; extend it by inserting extra midpoint $t$-points around the
  failure region.
- Pros: more robust for highly non-linear systems; builds directly on the
  existing $t$-grid logic; low risk.
- Cons: slower; more grid points to solve.
- Effort: low.
- Applicability: directly targets the $\phi_{HtM} = 0.2$ CD cell (which fails
  early, at $t = 0.70$) by densifying the solution path.

## Regularization

- Idea: add a small penalty term to the objective or residual to avoid
  singularities near bifurcation points.
- Pros: stabilizes the solver close to degenerate configurations.
- Cons: introduces bias into the solution; the penalty weight needs tuning.
- Effort: low-medium.
- Applicability: useful if the failure is near a genuine singularity (for
  example the complementarity boundary at high HtM).

## Alternative Solvers (trust-region / Levenberg-Marquardt)

- Idea: swap the root finder for a trust-region or Levenberg-Marquardt method
  with better ill-conditioning handling.
- Pros: robust to poorly conditioned Jacobians.
- Cons: must be adapted to the complementarity structure; rewrites part of the
  solver interface.
- Effort: medium.
- Applicability: a secondary option if Route B is unavailable.

## Adaptive Step Sizing

- Idea: adjust the $t$-grid spacing dynamically based on solve success or
  failure, inserting points where the solver struggles.
- Pros: steers around problematic regions of parameter space.
- Cons: too-large steps can miss important dynamics.
- Effort: low.
- Applicability: overlaps with enhanced homotopy; a good safety net.

## Preconditioning (variable rescaling)

- Idea: rescale variables to improve the numerical conditioning of the system.
- Pros: reduces sensitivity to floating-point error.
- Cons: requires domain-specific knowledge of the scales involved.
- Effort: medium.
- Applicability: a general robustness improvement that helps all cells.

## Better Initial Guesses and Neighbor Seeding

- Idea: seed each cell from a converged neighbor (already used) and add a
  damped Newton step or a line search.
- Pros: cheap; often sufficient on its own.
- Cons: limited if the solution path itself is unstable.
- Effort: low.
- Applicability: a quick win; combine with homotopy continuation.

# Recommended Sequencing

A pragmatic order that front-loads the cheapest, lowest-risk changes:

1. Quick wins: better neighbor seeding plus a denser or adaptive homotopy around
   the two failing $t$-ranges. Low effort and may already suffice.
2. Route B restoration on the host Mac: the highest-confidence principled fix.
3. If Route B is memory-infeasible: add regularization and preconditioning on
   top of Route A.

# Verification Plan

After any change, re-run the HtM sweep on the host and check:

- Both cells now reach $t = 1.0$ with the residual below tolerance.
- The reconstructed-vs-converged gap is quantified (especially the CD
  $\phi_{HtM} = 0.2$ underestimate).
- `VALIDATION.md` HtM table is updated to drop the "Reconstructed" flags.

```{.bash}
julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid()'
```

# Open Questions

- Can Route B run within the host Mac RAM budget (README cites 32+ GB)?
- Is the $\phi_{HtM} = 0.2$ CD failure a genuine economic singularity or a
  solver artifact? Diagnose via the residual norm and the Jacobian condition
  number at the last convergent $t$.
- Should the reconstructed values be excluded from Figure 4 until they are
  genuinely converged, to avoid presenting estimates as exact results?
