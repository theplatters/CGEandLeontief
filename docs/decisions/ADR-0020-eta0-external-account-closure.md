# ADR-0020 — The eta = 0 endpoint's external account: closure options

- **Status:** proposed (awaiting the user's decision; the audit that motivates
  it is settled, the choice is not)
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

Not taken yet — this record exists so the choice is made once, with the
measurements above in hand, and so the v6 scope is visible before any code is
written. The recommended sequence, for the user to confirm or reject:

1. **Adopt A now.** The v5 generation stays the citable one; the BF rows are
   read from the resource side and the gap; the `labor.BF` open gate records
   the non-identification. This is what the same-date corrections in
   `docs/DOCS_ASSESSMENT.md` (Version 7, section 4.2), the v5 flow table and
   the session log already do.
2. **Choose B or C before any v6 run.** B if the paper needs an *identified*
   BF external position and accepts BF = ALPHA at the aggregate level; C if
   the BF row must carry the immobile-benchmark economics, which is what the
   labour-closure narrative and `cbase2/review.md` section 1 ask for. B and C
   are not exclusive: B closes the account, C gives the allocation content,
   and a v6 could carry C with B's labour equation as the accounting closure.

A third possibility is recorded so it is not re-litigated: dropping the BF
row from the matrix and reporting ALPHA/GAMMA/DELTA only. It is not
recommended — BF is the paper's friction arm — but if C is out of scope, the
BF row currently duplicates ALPHA and the matrix loses nothing by saying so.

## Consequences

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
