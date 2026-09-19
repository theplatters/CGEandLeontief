# ADR-0021 --- Wage structure in the fixed-wage closure (GAMMA)

Status: proposed (2026-09-18)

Related: ADR-0002 (naming), ADR-0004 (run records), ADR-0005 (closure
plugins), ADR-0006 (preregistration), ADR-0014 (real wage, admissibility),
ADR-0017 (round-gain criterion), ADR-0019 (external account), ADR-0020
(sectoral wages at eta = 0, the vector-wage apparatus this builds on).

## Context

The fixed-wage closure pins one real wage as the numeraire, keeps the
allocation fully mobile, and drops the labour equation and the external
unknown. In the evaluation matrix its only free dimension is the financing
row, and that row is nearly degenerate: F2 and F3 are real-neutral by theorem,
and DELTA reproduces GAMMA to six digits. Of the five labour rows, GAMMA and
DELTA therefore contribute essentially one degree of freedom, and within GAMMA
the F2 and F3 columns are identical.

The flatness is a design property, not a modelling one. The closure's
parameter space already contains a wage vector; the design sets it to its
degenerate value (all sectors equal) and the kernel hard-codes that value as
the scalar `w = 1`. The question this ADR answers is whether the vector-wage
apparatus built for the eta = 0 endpoint (ADR-0020) can be reused to make the
wage structure an explicit, measurable dimension of the fixed-wage regime.

Measured evidence (full-71 A-bill, `GAMMA-F2` cell,
`experiments/probes/probe11_gamma_wage_structure.jl`, 2026-09-18):

- Pinning a vector reproduces the kernel cell at the degenerate point
  (`wbar = 1`: residual 8.9e-16, employment 1.00125981).
- A uniform rescale is a numeraire change: `wbar = 0.9` and `1.1` leave the
  allocation, employment and real income invariant to 8.9e-16 (the CPI scales
  exactly). The wage level is not an instrument.
- A relative change is a genuine shock: a ten percent wage push in the seven
  programme sectors moves employment by about half a percent (1.00686 against
  1.00126) and real income by about one percent (1.01186 against 1.00126);
  a single-sector tilt moves the allocation by about two percent in log terms.
- The external-account identity closes under every tilt: the price-weighted
  clearing residual stays at the solver-residual level (2.9e-16 at the
  degenerate pin, at most 7.0e-12 under a tilt), and that quantity equals the
  kernel canary's `diff` exactly along quantity perturbations (ratio 1.0000 at
  one and ten per mille). The proxy is therefore the gap; it was validated off
  equilibrium because both quantities are machine zero on it.
- Conditioning is unaffected: the Jacobian condition number stays between 6.32
  and 6.34 across all scenarios, and the admission criterion of ADR-0017
  involves only technology and demand shares, never the wage.

## Options

| Option | Mechanism | Cost | Assessment |
| --- | --- | --- | --- |
| A. Pinned wage vector as a parameter of the existing closure | `problem_fixed` / `_solve_fixed` take `wbar` (default: the scalar 1, so every current cell is bit-identical); the design pins `wbar` per cell | Dispatch-only change; the canary and `gdp_components` must accept the vector on their 2N branch; a new generation for provenance | Recommended: no new closure, the degenerate case is preserved by construction, and the wage structure becomes a design dimension |
| B. A separate closure id (for example `GAMMA-W`) | A new registry entry, its own formulation, tests and ADR | Registry and test surface grow for what is one parameter | Rejected: the closure is "a fixed real wage"; the scalar is its degenerate point, not a different closure |
| C. No kernel change; widen GAMMA through sectoral supply shocks only | A per-cell supply-shock spec in the design schema and `build_cell_model` (the kernel already takes a supply shock; the harness hard-codes a null one) | No kernel change, but a design-and-harness batch and a new generation | Viable, and it also restores the identifiability of `eta_s` and of the CES-versus-Leontief contrast; kept as the complement rather than a substitute |

## Decision (proposed)

Adopt option A. The fixed-wage closure takes an optional pinned wage vector,
defaulting to the current scalar pin, so that every existing GAMMA and DELTA
cell is unchanged by construction and the twelve non-BF cells of the current
generation reproduce their manifests. The wage structure is then a dimension of
the *design* (preregistered per cell), not a new closure, and the registry's
`labor.GAMMA` entry gains a sentence describing the vector formulation.

Two semantics are fixed by this ADR, because they are properties of the
closure rather than choices:

- The wage *level* is not an instrument. Any uniform rescale of the pin is a
  numeraire change, so scenarios must vary the structure. This holds because
  the model has no nominal anchor besides the wage; with a nominal government
  bill, nominal debt or money the invariance would fail, and the ADR should be
  revisited if such a block is added.
- The pinned vector is an exogenous scenario parameter, never an equilibrium
  object. Wage-structure cells are instrument responses, like the programme
  incidence or the shock size, and must be preregistered before the run that
  reports them.

## Consequences

- A `src/` change is required, so the current generation's provenance is
  superseded and the matrix must be re-run as a new generation before any
  wage-structure cell is cited (ADR-0004). The re-run is expected to reproduce
  the existing cells exactly, which is itself the promotion's test.
- `external_balance_canary` and `gdp_components` currently pin `w = 1` on their
  2N branch; both must take the vector, otherwise the acceptance gate does not
  run on the kernel's own arithmetic for a tilted pin. The probe's proxy
  measurement is evidence that the gate will pass, not a substitute for it.
- The reported `wage` metric stays a scalar (the wage-bill-weighted aggregate,
  as at eta = 0) with `wage_min` and `wage_max` carrying the structure.
- The F3 booking is priced at the equilibrium prices, so it moves with the wage
  structure; the identity gate is re-checked per scenario.
- Paper-facing language must keep the two limits: a wage-structure scenario
  does not identify a labour-supply elasticity, and a welfare difference across
  wage scenarios is not a statement about wage policy.

## Enforcement (proposed)

- The degenerate default reproduces the current cells: re-solving the twelve
  non-BF cells and the GAMMA/DELTA rows must match the committed manifests.
- A test asserting level invariance: the same cell under a uniformly rescaled
  pin returns the same allocation, employment and real income, with only the
  numeraire moving.
- A test asserting the identity under a tilt: `|S + T + M - (I+X) - (F + B_gov)|`
  at the solver-residual level for a tilted pin (measured worst 7.0e-12).
- A preregistered design carrying the wage-structure cells; the scenario set is
  the one open decision this ADR leaves to the operator (candidate set: a
  uniform rescale as the numeraire control, a programme-aligned push and a
  programme-aligned restraint, and a broad tilt as a robustness case).
