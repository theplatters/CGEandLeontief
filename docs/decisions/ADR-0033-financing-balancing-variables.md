# ADR-0033 — Financing formulations name the balancing variable

- **Status:** proposed
- **Date:** 2026-10-18
- **Supersedes:** — (amends the formulation text of ADR-0002 rows `F2` and
  `F3`; the closure set and the F1/F2/F3 taxonomy are unchanged)
- **Related:** ADR-0002 (closure taxonomy), ADR-0019 (endogenous external
  transfer `F`), ADR-0020 (the η = 0 endpoint carries sectoral wages),
  ADR-0032 (every labour row); `registry/closures.toml` `[financing.F2]`,
  `[financing.F3]`; `docs/closures/financing.md`; `revision/tables/table2_closure-conditions.tex`,
  `table3b_closure-matrix-3x5.tex`, `table3c_closure-matrix-3x5.tex`

## Context

ADR-0002 wrote the financing closures as budget identities:

- `F2`: `Σ p_i g_i = T(p)`, described as a "lump-sum / balanced-budget tax";
- `F3`: `Σ p_i g_i = F`, described as an "external balance".

Three things in these formulations do not match the implementation or the
accounting identities:

1. **The tax is not lump-sum.** The kernel (`tau_rate`, `OpenCGE`
   `src/closures/financing/financing.jl`) sets a proportional wage-income tax,
   `τ = Σ p g / (w·ΣL)`, and `E^h = (1−τ)·w·ΣL + F`. Its own comments call the
   base "lump-sum", which describes only the aggregate budget: with one
   household and no labour-leisure choice in the closure equations, a
   proportional wage-bill tax and a fixed amount give the same household budget.
   The two differ in incentives, which matters under the elastic labour supply
   (BETA).
2. **The tax also finances the baseline government purchases.** `tau_rate`
   finances `g^G` (the baseline government demand) in every closure, not only
   under F2. Under F2 the programme `g` is added to that base. The formulations
   and Tables 2, 3b and 3c stated F1 and F3 without any tax, which is wrong.
3. **`F` is not `Σ p g`.** ADR-0019 makes `F` the endogenous net external
   transfer that closes the external account
   `S + T_int + M − (I+X) = F + B_gov`. `F` is an unknown of the mobile
   solves (ALPHA, BF including η = 0, sectoral BETA) and is identically zero
   in the fixed-wage solve (GAMMA, DELTA). The registry records `F + B_gov` as
   the financing-neutral quantity at η = 0 (ADR-0020). The programme cost
   `B_gov = Σ p g` is booked externally under F3 only.

The paper tables carried these statements (Table 2, Table 3b, Table 3c), and
`docs/closures/financing.md`, the registry and the closure skill repeated them.

## Decision

1. **The tax is a proportional wage-income tax at an endogenous rate `τ`.**
   The word "lump-sum" is removed from the formulations. It may appear only to
   describe the aggregate household budget, which is the same as under a
   lump-sum tax in this model.
2. **`τ` finances the baseline government purchases in every closure:**
   `τ·w·ΣL = Σ p g^G`. Under F2 the programme is added to the base:
   `τ·w·ΣL = Σ p (g^G + g)`.
3. **`B_gov` is the primary symbol for the programme cost.** `B_gov = Σ p g`
   under F3, booked externally; `B_gov = 0` under F1 and F2.
4. **`F` is the net external transfer, fixed by the external account
   `S + T_int + M − (I+X) = F + B_gov`.** It is endogenous on the mobile labour
   rows and identically zero on the fixed-wage rows (GAMMA, DELTA). It is not
   defined as `Σ p g`.
5. **Household expenditure** is `E^h = (1−τ)·w·ΣL + F` in every closure, with
   `F = 0` on the fixed-wage rows.
6. **F1** is the preference reallocation `β̃_i ∝ β_i d_i` with budget-neutral
   household expenditure. It has no additive demand and `B_gov = 0`.

The change is to formulation text and to the tables and documents that quote
it. The kernel and the run manifests are unchanged, so no pin bump and no new
generation are needed (the implementation already does this).

Alternatives considered:

- *Keep `T(p)` and the lump-sum label.* Rejected: it misstates the tax base
  and the incentive channel, and omits the baseline government financing.
- *Rename to a wage-tax label without changing the formulas.* Rejected: the
  formulas were also wrong (F1 and F3 without `τ`, `F` defined as `Σ p g`).

## Consequences

- Tables 2, 3b and 3c state the same conditions (the revised text is in the
  files above). Table 3c's caption defines `E^h`, `τ`, `F` and `B_gov` once.
- `registry/closures.toml` `[financing.F2]` and `[financing.F3]` carry the new
  `formulation` strings; `docs/status.md` shows them and must be regenerated
  (`scripts/status.jl`).
- `docs/closures/financing.md` and the closure skill (`.claude/` and
  `.opencode/` copies, byte-identical) are corrected.
- Not yet corrected, and to be done in follow-up commits: `docs/DOCS_ASSESSMENT.md`
  (lines ~291, 301, 343, 710–711; this file needs a version bump and revision
  colouring under the `md-style-corrections` scheme), the paper's model section
  `revision/sections/(B)model-closures-proofs.tex` (lines ~276–286, which still
  describes a lump-sum tax `T(p)`), `revision/sections/(F)closures-mapping.tex`
  (lines ~143–147), and the kernel docstring in `OpenCGE`
  (`src/closures/financing/financing.jl`, the phrase "lump-sum tax base"), which
  is a kernel-side wording change and goes through `../opencge`.
- ADR-0002 stays `accepted`. Its F2 and F3 formulation text is amended by this
  record; its closure set is unchanged.

## Verification

Kernel paths are those of the pinned OpenCGE package (`v0.2.1`, the copy
`Manifest.toml` resolves):

- `tau_rate` (`src/closures/financing/financing.jl`, lines 120–123): the tax
  base is `gov_demand` for every financing, `gov_demand + g` under F2.
- `household_expenditure` (lines 140–145): `E = (1 − τ)·w·ΣL + F`.
- `external_balance` (lines 177–182): zero except under F3.
- `solve` (`src/core/equilibrium.jl`): the fixed-wage path passes
  `external_transfer = 0.0` (line 1115, 1130); the mobile path takes
  `F = x[2N+2]` (η = 1) or `x[3N+1]` (η = 0 and sectoral BETA, line 1257).
- The external canary books `B_gov = dot(p, additive)` under `ExternalDebt`
  only (line 805) and reports `S + T + M − IX − (F + B_gov)` (line 810).

Cell checks on the real table (full-71 calibration, programme `g = 0.001` per
sector, `θ = ϵ = 0.5`, `σ = 0.9`, labour elasticity `η` as in the row; the
script is kept outside the repository):

| Labour | Financing | F | τ vs formula | E^h vs (1−τ)wΣL + F | B_gov | External canary |
| --- | --- | --- | --- | --- | --- | --- |
| ALPHA (η = 1) | F1 | 0 | exact | exact | 0 | 2e-16 |
| ALPHA (η = 1) | F2 | −7.4e-3 | exact | exact | 0 | 9e-17 |
| ALPHA (η = 1) | F3 | −7.8e-2 | exact | exact | 0.071 | 2e-16 |
| BF (η = 0) | F1 | 0 | exact | exact | 0 | 1e-16 |
| GAMMA (fixed) | F1 | 0 | exact | exact | 0 | not run |
| GAMMA (fixed) | F2 | 0 | exact | exact | 0 | not run |
| GAMMA (fixed) | F3 | 0 | exact | exact | 0.071 | not run |

Not verified in this run: BF (η = 0) under F2 and F3. With this placeholder
programme the solve did not converge (`Stalled` and `MaxIters`, residual
5e-4); the registry records these cells as executed under the production
programme bundle, and this ADR does not re-verify them. The GAMMA canary was
not re-run; the registry reports GAMMA-F3 closing to about 1e-12.

Two further points from the run:

- `τ` is an equilibrium variable, not a constant. It equals the baseline 0.2141
  in ALPHA-F3 but 0.1951 in GAMMA-F3, because the fixed-wage employment moves
  the tax base `wΣL`.
- The earlier draft of this ADR, and of Tables 2, 3b and 3c, called `F`
  endogenous in every closure. The fixed-wage solve sets it to zero, and this
  section corrects that.

## Enforcement

- `scripts/check_repo.jl` must report 0 violations before and after the commit
  that lands this ADR (AGENTS.md, Runs).
- A later change to the `τ` base, the `F` definition or the `B_gov` booking is
  a new ADR, not an edit to this one.
