# ADR-0020 — The eta = 0 endpoint's external account: closure options

- **Status:** accepted (user instruction, 2026-09-18): **option C** — sectoral
  wages at the eta = 0 endpoint — with option B retained as a future venue and
  option A superseded by this choice. Specification and probe evidence below;
  the audit in Context is unchanged.
- **Date:** 2026-09-18
- **Supersedes:** —
- **Related:** ADR-0002, ADR-0004, ADR-0006, ADR-0010, ADR-0019;
  `src/core/equilibrium.jl` (`problem`, the eta = 0 pin), `src/core/diagnostics.jl`
  (`external_balance_canary`, `gdp_components`), `tests/test_external_closure.jl`,
  `experiments/probes/probe5_bf_f3_external_position.jl`,
  `experiments/probes/probe5b_bf_f3_followups.jl`,
  `experiments/probes/probe5c_bf_pin_equivalence.jl`;
  `registry/closures.toml` entries `labor.BF` and `financing.F3`;
  `paper/tables/matrix_5x3_v5_flows.md` (the correction of the same date);
  `cbase2/review.md` section 1 (the sector-specific-wage point)

## Context

ADR-0019 gives every regime all N goods-market clearings and one scalar
external unknown where it is needed. At the BF eta = 0 endpoint it does not
introduce an unknown: the sectoral allocation is frozen at
`data.labor_share`, so the labour-supply residual is identically zero
(`sum L_i = sum labor_share = Lbar` for every trial point) and the closure
replaces that row with the explicit pin `F = 0` (ADR-0019 decision D2). The
pin is therefore a *substitute for the labour-market equation*, not a
statement about the external account. That has three consequences, all
measured on the executed v5 generation on 2026-09-18:

1. **The reported position is the booking.** `external_position = F + B_gov`
   collapses to `B_gov = sum(p g)` under the pin. It equals the programme
   cost, 1.3310 percent of GDP (40 300 EUR m), and is identical to the last
   digit in BF-F3, GAMMA-F3 and DELTA-F3 — a quantity that no labour closure
   can move is not a labour-market result. GAMMA/DELTA-F3 are nonetheless
   genuine positions: their account closes to ~1e-12 (employment rises to
   1.017796 and the resources arrive from abroad), while BF-F3's does not.
2. **The account is open, and the gap is the pin's labour residual.** The
   resource side `S + T_int + M - (I+X)` is +0.5374 percent of GDP at BF-F3
   against the booked +1.3310 percent, the gap -0.7936 percent; the identity
   gap is exactly `-w * (sum L^cm - sum L)` with `sum L^cm = 1.0079361`
   against the frozen bar 1 (verified at all 15 v5 cells; at the twelve
   eta = 1 cells the two quantities agree to machine precision). So the
   frozen allocation is not what creates the gap: the gap closes as soon as
   the pin is set to the labour-clearing value.
3. **The position is not identified at eta = 0.** Replacing the pin by
   `F = c` leaves every `c` an exact root (residual 2.2e-15) and moves the
   reported entry one-for-one: `c` in {-0.02, -0.01, 0, +0.01, +0.02} gives
   -0.669 / +0.331 / +1.331 / +2.331 / +3.331 percent of GDP, with the
   resource side moving too (-0.355 / +0.091 / +0.537 / +0.983 / +1.429).

Two further measurements bound the options. Solving the eta = 0 system with
the pin set to ALPHA's solved `F` reproduces the ALPHA cell to 2.8e-17 in the
canonical vector, in household expenditure and in the consumption block; and
replacing the pin by the labour equation `sum L^cm = Lbar` closes the account
to ~4e-16 and returns BF-F1/F2/F3 = ALPHA-F1/F2/F3 exactly, with financing
neutrality then holding at eta = 0 as well (`F_F3 = F_F2 - B_gov` to 4e-17).
The reason is structural: in this demand-only, one-factor, CRS design the
frozen sectoral allocation enters nothing aggregate — only `sum L_i`, which
equals the bar either way. The entire BF/ALPHA difference in v5, including
the BF-F3 welfare entry of 0.000 percent (which becomes -1.8226 percent at
the labour-clearing `F`), is carried by the pin.

This converges with two already-recorded open gates: the `labor.BF` gate
"eta = 0 keeps one common wage with fixed quantities; BF's immobile case has
sector-specific wages" (`cbase2/review.md` section 1) and the session log's
"the eta = 0 sector-specific-wage gap". A closure that prices the frozen
allocation would give the labour-closure axis a real eta = 0 arm; the current
pin does not.

## Options

| Option | Mechanism | Measured effect | Cost and risk |
| --- | --- | --- | --- |
| A. Keep the pin; report the position as not identified (documentation only) | `F = 0` pin stays (ADR-0019 D2); the BF rows report the resource-side imbalance and the gap; no external position is claimed at eta = 0 | Account open by -0.7936 percent of GDP at BF-F3; the entry stays a design constant | No re-run; v5 remains citable. The BF row then carries no external position at all, and the "immobile benchmark" is a normalisation of ALPHA rather than a distinct economy |
| B. Labour equation instead of the pin at eta = 0 | `out[2N+1] = sum L^cm(p, y, w) - Lbar`, `F` free | Account closes to ~4e-16; `F` identified; BF-F1/F2/F3 = ALPHA-F1/F2/F3 exactly; financing neutrality extends to eta = 0 | `src/` change, so a new generation (ADR-0004). BF collapses onto ALPHA in every aggregate, so the row adds no aggregate information in a demand-only design — it must be paired with a design in which the allocation bites (supply-side shock, second factor) or the row is redundant |
| C. Sector-specific wages at eta = 0 (the true BF immobile benchmark) | No common wage: `w_i` from sector `i`'s own marginal-product condition at the frozen `L_i`, with the zero-profit block evaluated per sector | Not measured; needs implementation. Gives the allocation a price channel, so a sectorally concentrated shock can create scarcity and move prices | Largest scope: new kernel feature (a wage vector, the zero-profit and labour blocks rewritten for eta = 0), new tests, new generation. Changes every BF cell and the current BF/ALPHA nesting claim; GAMMA/DELTA/BETA are unaffected if the change is confined to eta = 0 |

## Decision

**Option C is accepted: sectoral wages at the eta = 0 endpoint.** The
specification below was derived and then validated on the full-71 A-bill
calibration before any `src/` change (probe evidence in the same section);
promotion into the kernel still requires the ADR-0006 workflow (tests, a
`matrix_5x3_v6` generation, registry rows) because a `src/` change invalidates
the provenance of every existing run.

### Specification

The eta = 0 endpoint keeps the frozen sectoral allocation `L_i = labor_share_i`
and replaces the single economy-wide wage by a wage per sector, set by that
sector's own marginal product at the frozen allocation:

- Unknowns: `X = [p(1:N); y(1:N); w(1:N); F]` — 3N + 1.
- Equations: (1) zero-profit per sector, `p_i = cost_i(p, w_i)` (N); (2) the
  sectoral first-order condition at the frozen allocation,
  `log L^cm_i(p_i, y_i, w_i) = log labor_share_i` (N); (3) all-N clearing, with
  household wage income `sum_i w_i L_i` and `E = (1 - tau) sum_i w_i L_i + F`
  (N); (4) the CPI numeraire (1).

`F` is kept, and this is the one substantive addition the counting forces: the
block {zero-profit, FOC, clearing} is homogeneous of degree 1 in `(p, w, F)`,
so its 3N equations determine 3N - 1 effective unknowns. Replacing the mobile
system's single aggregate labour equation by N sectoral conditions adds N - 1
equations, so the block has one equation more than it has directions to pin;
the demand block needs one free scalar to be consistent with the supply side.
`F` is that scalar, entering `E` after tax exactly as in ADR-0019. The
alternative (drop one clearing equation) is the ADR-0010 shortcut that
ADR-0019 retired, so it is not available.

Two properties follow, both measured:

- **The external account closes at eta = 0.** At a FOC solution the frozen
  allocation is cost-minimizing at these wages, so the identity gap
  `-(sum_i w_i L^cm_i - sum_i w_i L_i)` vanishes. The factor-market gap that
  the `F = 0` pin left (up to -0.7936 percent of GDP) is gone by construction.
- **Financing neutrality extends to eta = 0.** With `F` free, the same algebra
  as ADR-0019 decision item 7 applies at the immobile endpoint: F2 and F3 have
  identical real allocations and `F_F3 = F_F2 - B_gov`, hence identical booked
  positions. The pin had destroyed this; the sectoral-wage closure restores it
  in both regimes.

### Probe evidence (2026-09-18, full-71 A-bill, `experiments/probes/probe7_sectoral_wages_eta0.jl`)

| Cell | max abs residual | identity gap | `F` | `B_gov` | booked `F + B_gov` |
| --- | ---: | ---: | ---: | ---: | ---: |
| BF-F1 | 1.4e-13 | +1.9e-14 | -0.00581497 | 0.00000000 | -0.00581497 |
| BF-F2 | 4.2e-13 | -4.2e-13 | -0.00852841 | 0.00000000 | -0.00852841 |
| BF-F3 | 1.7e-12 | -5.0e-13 | -0.02337469 | +0.01484628 | -0.00852841 |

Also measured: the FOC gap `max|log L^cm - log L|` = 1.7e-12 at BF-F3;
`F_F3 = F_F2 - B_gov` with BF-F2 and BF-F3 identical in allocation and booked
position (neutrality at eta = 0); `F` identified with an F-column norm of 0.163
in the finite-difference Jacobian; multi-start invariance at 2.1e-9 in the
unknown vector (dF ~ 4e-12). Two caveats to carry into promotion: the Jacobian
is stiffer than the mobile all-N system (condition number 9.4e7,
`sigma_min/sigma_max` = 1.1e-8, against 52.8 for the ADR-0019 mobile system), so
the acceptance gate and the polish target need to be set from measurement; and
`B_gov` is priced at the eta = 0 equilibrium prices, which are not one (the
sectoral wages move them), so the F3 booking reads +1.4846 percent of GDP here
against the baseline +1.3310 percent. Sectoral wages at the baseline-type cell
span 0.970 to 1.672 (ratio 1.7) around the CPI-pinned level, and the wage bill
`sum_i w_i L_i` = 1.00702.

### Status of the other options

- **A (keep the pin, report the position as not identified)** is superseded by
  this choice: it was the documentation-only response to the audit, and C
  removes the problem it documented.
- **B (labour equation instead of the pin, `F` free, common wage)** is retained
  as a future venue, not adopted: it closes the account and identifies `F`, but
  it collapses the BF row onto ALPHA in every aggregate (measured to 2.8e-17),
  so it buys accounting coherence by giving up the immobile benchmark as a
  distinct economy. It remains the cheapest closure if a generation needs a
  closed account at eta = 0 without new kernel machinery; a future ADR can
  supersede this decision to adopt it.
- **Dropping the BF row** is not adopted: C gives the row its own economics.

### Kernel changes required for promotion (not yet made)

`_cost_minimizing_labor` must accept a wage vector (it hard-codes `log(w)` for
a scalar); `_mobile_market_demand` needs the sectoral-wage variant (frozen
`L_i`, income `sum_i w_i L_i`); `problem` needs the eta = 0 sectoral branch
(3N + 1 unknowns, the four blocks above); `external_balance_canary` and
`gdp_components` must accept the wage vector; `Solution.wages` already stores a
vector, so the reporting layer needs no new field. `tests/test_external_closure.jl`
currently asserts the pin and the measured pin-gaps at the BF cells and must be
rewritten for the new contract (account closed to <= 1e-12 at eta = 0, `F`
identified, neutrality at eta = 0).

## Consequences

Operative choice: **C** (decision above). The A and B lines below are kept as
the record of the alternatives considered.

- Under A, no kernel, test, run or registry status changes; only the
  documentation of the reading changes, and `registry/closures.toml`
  (`labor.BF`, `financing.F3`) records the non-identification as an open
  gate. `matrix_5x3_v5` remains the citable generation.
- Under B or C, `src/` changes invalidate the provenance of every existing
  run (AGENTS.md), so the matrix is re-minted as `matrix_5x3_v6` with a new
  preregistration (ADR-0006) and the BF rows are re-derived. Under B the
  BF-F1/F2/F3 manifest values change (external transfer, position, gap) and
  the BF testset in `tests/test_external_closure.jl` must be rewritten: it
  currently asserts the `F = 0` pin and the measured gaps
  (+4.3438e-04 / -5.6182e-04 / -7.9361e-03), which B removes by construction.
- Under C the BF/ALPHA nesting (ALPHA = BF at eta = 1) survives only if the
  wage-vector change is confined to eta = 0; the BF row then becomes the only
  arm in which the allocation has a price channel, which is the point.
- In every case the paper-facing tables must stop printing a BF external
  position: the v5 flow table already carries the resource side and the gap,
  and the corrected Version 7 text states the non-identification.

## Enforcement

- A (adopted now): the `F = 0` pin stays asserted
  (`tests/test_external_closure.jl`, the BF testset) and the eta = 0 identity
  gap stays reported, not gated; the `labor.BF` and `financing.F3` open gates
  in `registry/closures.toml` name the non-identification and this ADR.
- B: a new eta = 0 gate asserting `|S + T_int + M - (I+X) - (F + B_gov)| <= 1e-12`
  at the BF cells (replacing the "reported, not gated" clause), the absence of
  the pin, and financing neutrality at eta = 0
  (`F_F3 = F_F2 - B_gov`, identical positions).
- C: a fixture in which the frozen allocation bites — a sectorally
  concentrated shock that moves prices under sector-specific wages while the
  common-wage fixture does not — plus the existing aggregate-labour and
  identity gates.

Either of B or C is a new ADR superseding this one's Decision section; the
measurements in Context stand as recorded evidence for both.

## Implementation (promoted 2026-09-18)

Option C is in the kernel and executed as the `matrix_5x3_v6` generation. The
change is dispatch-only, so every η = 1 and fixed-wage path is bit-identical:

- `_cost_minimizing_labor` evaluates `log.(w)`, so it accepts a wage vector and
  reduces bit-identically to the scalar case; `_wage_bill(w::Real, L) = w·ΣL`
  keeps the scalar arithmetic exact, while `w::AbstractVector` sums `w_i·L_i`.
- `problem_sectoral` is the η = 0 system (3N+1 unknowns `[p; y; w(1:N); F]`: N
  zero-profit, N sectoral FOC `log L^cm_i = log L̄_i`, N clearing, one CPI
  numeraire). `problem` now rejects η = 0 (the pin branch is gone);
  `equilibrium_residuals`, `market_clearing_residuals` and
  `external_balance_canary` accept the 3N+1 vector and expand a legacy 2N+2
  vector by replicating the scalar wage.
- `solve` dispatches `problem_sectoral` at η = 0 with a tighter primary
  tolerance (1e-8) and a longer polish ladder (6 steps, target 1e-13), because
  the system is stiff (Jacobian condition ~9.4e7 against ~52.8 mobile); the
  η = 1 settings are unchanged.
- `gdp_components` uses the sectoral wage vector at η = 0; the harness reports
  the wage-bill-weighted average as the scalar `wage` metric and the dispersion
  as `wage_min` / `wage_max`; the η = 0 "third gate" is the sectoral
  labour-market gap (`sectoral_labor_gap`, exported), and
  `assert_external_account` now gates the identity in every regime.

Enforcement as realised: the pin is absent, financing neutrality holds at
η = 0 (`F_F3 = F_F2 − B_gov` to 3.5e-14), and the identity is gated. The gate
value deviates from the 1e-12 written in the Enforcement section above: the
measured floor of the stiff η = 0 system is 1.04e-11 (BF-F1, against 1.8e-13 at
BF-F2/F3), so the harness gate is 1e-9 and the test assertion 1e-10 — set from
the measured floor, not tuned to pass a failing cell.

Results (15/15 cells executed, commit `6db1da5`):

| Cell | `F` | `B_gov` | Net ext. pos. | Gap | Consumption rel. | Wage vector |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| BF-F1 | -0.005815 | 0 | -0.005815 | -1.0e-11 | +0.096 % | 0.9734 … 1.5249 |
| BF-F2 | -0.008528 | 0 | -0.008528 | -1.8e-13 | -1.980 % | 0.9702 … 1.6719 |
| BF-F3 | -0.023375 | +0.014846 | -0.008528 | -1.9e-13 | -1.980 % | 0.9702 … 1.6719 |

v5 for comparison: net external position 0 / 0 / +1.331 % of GDP, gap
+4.3e-04 / -5.6e-04 / -7.9e-03, consumption rel. +0.043 % / -1.694 % /
0.000 %.

Verification: the twelve non-BF cells reproduce v5 (ALPHA/BETA to 6.7e-16,
GAMMA/DELTA to 1.6e-11 — the fixed-wage warm-start sensitivity (well-conditioned Jacobian, cond ≈ 6.3; the metric layer amplifies the residual level), which
appears identically on the pristine kernel, so it is not an effect of this
change). Evidence: `experiments/probes/probe7_sectoral_wages_eta0.jl` (closure),
`probe8_promotion_verification.jl` (kernel and harness),
`probe9_nonbf_reproduction.jl` (baseline comparison),
`probe10_repin_eta0_fixtures.jl` (fixture goldens); runs in
`runs/matrix_5x3-v6-*`.
