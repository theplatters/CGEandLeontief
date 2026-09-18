# ADR-0019 — All-N goods-market clearing with an explicit external account

- **Status:** accepted (user instruction, 2026-09-18)
- **Date:** 2026-09-18
- **Corrects:** ADR-0010's mobile N−1 formulation decision (the N−1 clearing
  decision only). ADR-0010 is NOT superseded as a whole: its kept BF endpoints
  η ∈ {0, 1} and the retired interpolation/allocation wedge stand.
- **Related:** ADR-0004, ADR-0005, ADR-0010, ADR-0012, ADR-0013, ADR-0016,
  ADR-0018; ROADMAP.md Phase-4 omitted-equation-invariance/goods-market-clearing
  gates; `src/core/equilibrium.jl`, `src/closures/financing/financing.jl`,
  `src/core/diagnostics.jl`, `tests/test_external_closure.jl`,
  `experiments/run.jl`

## Context

The kernel enforces N−1 goods-market clearings plus the CPI = 1 numeraire in
the mobile-labour system and exposes the omitted N-th residual as the
"external-account canary" (ADR-0010). ADR-0010 rejected the `revisefinal` all-N
mobile form as over-determined for open economies: N clearings without the CPI
numeraire leave the residual homogeneous of degree 1 in `(p, w)`, so the
equation count exceeds the effective unknown count by one. Measured at the
time, the all-N solve stalled at a least-squares floor (≈ 4.4e-4 on the real
calibration; 9.0e-3 (F3) / 4.7e-4 (F2) on the v3 contract fixture).

That diagnosis is now understood to have missed the variable, not the
equation: the all-N system lacked the external account as an unknown. With one
scalar external transfer added, the all-N system is square and has an exact
root. Probe evidence measured by the orchestrator on the full-71 A-bill
calibration: the baseline all-N residual at `(p, y, w, F) = (1, λ, 1, 0)` is
4.4e-16; the finite-difference Jacobian has condition number 52.8 with F-column
norm 0.163 (F is identified, not a null direction). ALPHA-F1/F2/F3 solve at
max|resid| ≈ 2.2e-15 with identity gap ≈ 2e-16; measured
`F_F2 = −1.014e-3`, `F_F3 = −1.4323e-2`, `B_gov(F3) = +1.3310e-2`, so
`F_F3 + B_gov = F_F2` with quantity agreement to 1.4e-16. The η = 0 BF cells
solve with the pin at ~2e-15, and the closed 2-sector tiny fixture
(s = m = X = I = 0) solves exactly with F = 0.

## Decision

1. **All N clearing equations are enforced in every regime.** No market is
   omitted; there is no residual market.
2. **The mobile η = 1 system gains one scalar unknown F** (net external
   transfer), entering household expenditure after tax:
   `E = (1−τ)·w·ΣL + F`. The solved system has 2N+2 unknowns `[p; y; w; F]`
   with N zero-profit + N clearing + labour-supply + CPI = 1 equations at
   η = 1. At the η = 0 (BF) endpoint the labour-supply row is identically zero
   and is replaced by the explicit pin row `F = 0`.
3. **The fixed-wage (GAMMA/DELTA) and η = 0 systems clear all N markets
   without F.** Their quantity margins (employment / the F = 0 pin) close the
   account; no external unknown is added there.
4. **External programme financing is booked as `B_gov = Σ_i p_i g_i` under the
   ExternalDebt (F3) closure; `B_gov = 0` otherwise** (F1 compositional, F2
   lump-sum programme tax).
5. **Exact identity at every all-N solution:**
   `S + T_int + M − (I+X) = F + B_gov`, where `S = s·E` with E including F,
   M is the full import content (final-import margin plus the
   intermediate-import leak `M_int`, row 74), and `T_int` (row 75) stays in T
   as before (ADR-0013).
6. **Baseline:** `F = 0` reproduces the calibration to machine precision; no
   recalibration.
7. **FINANCING NEUTRALITY (theorem of this closure).** Because F is an
   unrestricted after-tax household transfer, F2 and F3 have identical real
   equilibria `(p, y, w)`, with `F_F3 = F_F2 − B_gov` and identical net
   external position `F + B_gov`. *Algebra:* write baseline after-tax labour
   income `Y_d = (1−τ)·w·ΣL` and programme cost `G_prog = Σ_i p_i g_i`. Under
   F2 the programme is taxed domestically (`B_gov = 0`) so
   `E_2 = Y_d − G_prog + F_2`; under F3 the household is untaxed and the
   programme is booked externally so `E_3 = Y_d + F_3` with `B_gov = G_prog`.
   The N clearing equations see only total final demand `C(E) + g + …`; hence
   `E_2 = E_3` implies identical clearing residuals, i.e.
   `F_3 = F_2 − G_prog = F_2 − B_gov`, and all remaining real equations
   (zero-profit, labour, CPI) coincide. So `(p, y, w)` coincide and
   `F + B_gov` is invariant. F1 (compositional) remains distinct. The v1–v4
   F2/F3 differences were artifacts of the omitted-market shortcut, not
   economics.

## Consequences

- The canary `diff` (`external_balance_canary`) becomes the identity gap
  `S + T_int + M − (I+X) − (F + B_gov)`, evaluated with F threaded into S. It
  is ≈ 0 at every η = 1 solution (mobile ALPHA/BETA and fixed-wage
  GAMMA/DELTA); at the BF η = 0 endpoint the all-N clearings hold but the
  canary retains the documented fixed-allocation/factor-market gap (measured
  4.3e-4 / −5.6e-4 / −7.9e-3 on full-71 F1/F2/F3 — zero-profit prices the
  cost-minimizing labour demand, not the frozen baseline allocation), and the
  gap is reported, not gated.
- `gdp_wedge` (ADR-0018) is exactly `−canary.diff`; at η = 1 solutions it is
  ≈ 0, so the income- and expenditure-side real GDP coincide (`F` is booked
  inside `ΣV` through household expenditure, so no deflator adjustment is
  needed). The booked external position is reported separately as
  `financing = F + B_gov` (canary) and `external_financing` (gdp_components);
  at BF η = 0 cells `gdp_wedge` is the factor-market gap instead.
- `Solution` gains `external_transfer`; the solved mobile vector is 2N+2 with
  the F component last.
- Because `src/` changes invalidate run provenance (AGENTS.md), the executed
  v4 generation (ADR-0018, already run) is superseded for citation by a new
  `matrix_5x3_v5` generation carrying this closure. The plan had originally
  intended to bundle this closure with the v4 generation, but v4 had already
  executed; v4's runs stay visible as history (ADR-0004).
- ADR-0010's claim that the all-N `revisefinal` form is over-determined is
  corrected: it lacked the external variable, not an equation. The remainder
  of ADR-0010 (η ∈ {0, 1} endpoints, retired wedge) is unaffected.

## Enforcement

- New `tests/test_external_closure.jl` contract: all-N residual ≤ 1e-10 at
  every solved cell; identity gap `|S + T_int + M − (I+X) − F − B_gov|`
  ≤ 1e-12 at every η = 1 solution (mobile and fixed-wage); the BF η = 0
  factor-market gap is asserted as a finite, documented, non-gated quantity;
  baseline reproduction at F = 0 to machine precision; F identification
  (non-degenerate F column) and multi-start invariance; the η = 0 `F = 0`
  pin; closed-fixture behaviour (exact solve with F = 0);
  financing-neutrality equality `F_F3 = F_F2 − B_gov` with identical `(p, y,
  w)` and identical `F + B_gov`.
- Updated `tests/test_calibration.jl` (baseline/F = 0 reproduction),
  `tests/test_kernel_regression.jl` (goldens re-pinned to the 2N+2 all-N
  solutions), and `tests/test_gdp_measurement.jl` (income/expenditure
  coincidence with the booked `F + B_gov` wedge).
- Experiments gate in `experiments/run.jl` (`assert_external_account`): every
  cell asserts the all-N clearing residual; η = 1 cells additionally assert
  the identity gap; the manifest records `external_transfer`,
  `programme_financing`, and `external_position = F + B_gov`.
