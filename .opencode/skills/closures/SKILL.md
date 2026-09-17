---
name: closures
description: Use when adding, changing, promoting, or retiring a labour-market or financing closure in the BeyondHulten repo (BF, ALPHA, BETA, GAMMA, DELTA, future ZETA; F1, F2, F3) — registry entries in registry/closures.toml, status transitions with evidence, implementation in the single root src/ kernel, tests, and the Phase 2 promotion from the frozen cbase2 snapshot. Covers the eta vs eta_s symbol rule and the DELTA corner rule from ADR-0002.
---

# Closures: adding, changing, promoting

## Stable ids and formulations (ADR-0002)

Full taxonomy: `docs/decisions/ADR-0002-closure-taxonomy.md`. Never reuse an id for a different formulation; adding an id is an ADR-level change.

| Axis | ID | Formulation | Rule |
| --- | --- | --- | --- |
| Labour | `BF` | `L_i = L_fixed_i^(1-η) · L_costmin_i^η`, `Σ L_i = L̄`, `η ∈ [0,1]` | `η` is the BF reallocation friction, not a supply elasticity |
| Labour | `ALPHA` | `Σ L_i = L̄`; one economy-wide flexible wage; full cost-minimizing allocation | the `η=1` limit of BF **by construction** — not an independent mechanism |
| Labour | `BETA` | `Σ L_i = L̄ · ((w/P)/(w₀/P₀))^η_s` | `η_s` is the supply elasticity; never write it as `η` |
| Labour | `GAMMA` | `w/P = w̄` (=1); employment endogenous and uncapped | two-sided peg, not the one-sided wage floor |
| Labour | `DELTA` | `GAMMA` + Leontief limit (`θ, ϵ, σ → 0⁺`) at `η = 1` | a **corner** reached as GAMMA + Leontief technology; the endpoint row of the 5×3 matrix, not an independent mechanism |
| Labour | `ZETA` | `0 ≤ L̄ − L ⊥ w/P − ω̄ ≥ 0` | future (`idea`); needs an MPEC/smoothing formulation — plain square solvers cannot express it |
| Financing | `F1` | `β̃ᵢ = βᵢdᵢ / Σⱼβⱼdⱼ`; `Σ pᵢcᵢʰ = Eʰ`; budget-neutral | demand composition, never an autonomous investment multiplier |
| Financing | `F2` | `Σ pᵢgᵢ = T(p)`; lump-sum/balanced-budget tax | needs a balanced-budget counterparty |
| Financing | `F3` | `Σ pᵢgᵢ = F`; external balance | accounting open-economy closure, not a monetary mechanism |

Symbol discipline: `η` is always the BF reallocation parameter, `η_s` is always the BETA supply elasticity. They never appear unqualified as "the elasticity".

## One kernel rule (ADR-0001)

- Implement in the root `src/` only. Current homes: `src/mobile_labor.jl` (BF, ALPHA, GAMMA) and `src/interface.jl` (closure types). `docs/log/2026-09.md` names `src/closures/` as the Phase 2 promotion target — planned, not created yet.
- `cbase2/` is a frozen snapshot: read `cbase2/src/closures.jl` and `cbase2/src/financing.jl` freely and copy the economics from them, but never edit `cbase2/` and never build a parallel kernel.
- Closures are plug-in options registered centrally; do not duplicate the CES core.

## Add or change an entry in registry/closures.toml

Schema: `registry/README.md` (`[labor.<ID>]` / `[financing.<ID>]`). Required fields:

- `id`, `name`, `formulation` (exact string), `interpretation`, `status`
- `symbols` (grep-able Julia types/functions), `files`, `references`, `tests`, `evidence` (repo-relative paths)
- `dead_ends`, `adrs`, `open_gates`, `notes`

Checklist:

1. If the id or formulation is new, add/extend an ADR first; never redefine an existing id.
2. Implement in root `src/`; add `tests/test_<feature>.jl` and include it from `tests/runtests.jl`. Cover numerical results with `isapprox`, not exact equality.
3. Fill the registry fields with real paths; keep unresolved issues in `open_gates`, not prose.
4. Run `julia --project=. -e 'using Pkg; Pkg.test()'`.
5. Update evidence and status, then `julia --project=. scripts/status.jl` → 0 warnings → `julia --project=. scripts/status.jl --check`.
6. Commit as `closure(<id>): …` and log the session in `docs/log/YYYY-MM.md`.

`scripts/status.jl` warns on missing `files`/`tests`/`evidence` paths, `symbols` not found in the listed files, unknown scenario ids referring to the closure, and missing `DE-`/`ADR-` records.

## Status transitions (evidence gates)

`idea → spec → implemented → tested → validated`, plus:

- `rejected` — must link a `docs/dead-ends/DE-*` record.
- `superseded` — must name the replacing closure id.

Evidence per status (ADR-0003, `registry/README.md`): code in a listed file that runs (`implemented`); listed tests cover the contract and the suite passes (`tested`); the Phase 4 gates of `ROADMAP.md` pass and a recorded run manifest is linked (`validated`). Never promote a status without the corresponding evidence.

## Promotion from frozen cbase2 (Phase 2, planned)

- Sources: `cbase2/src/closures.jl` (BETA, DELTA) and `cbase2/src/financing.jl` (F1–F3). The `open_gates` for each entry list what blocks promotion.
- Port the economics into the root kernel, add contract tests, then record the promotion by updating the closure entry (`files`, `tests`, `status`, `open_gates`, `notes`) and the lab log.
- Never edit `cbase2/`; its freeze commit, tree hash, kernel diffs, and open items are in `registry/freeze.toml` (`[frozen.cbase2]`).
- Do not claim `tested`/`validated` for promoted code until the root tests actually cover it and pass.
