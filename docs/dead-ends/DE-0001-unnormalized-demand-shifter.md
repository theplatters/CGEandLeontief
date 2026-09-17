# DE-0001 — Unnormalized household demand shifter

- **Status:** dead end — do not retry without meeting the revival conditions
- **Recorded:** 2026-09-17
- **Origin:** rejected mobile-labour system (`src/mobile_labor.jl` pre-fix), 2026-09
- **Related:** ADR-0002, DE-0002, closures `BF` and financing `F1`
- **Evidence:** `src/mobile_labor.jl:148` (pre-fix line), `roadmaps/vertdict.md`,
  `ROADMAP.md` §3; diagnostic `Σ_i β_i d_i p_i^{1−σ} ≈ 1.25485`

## What was tried

Represent a demand shock by multiplying the household's CES weights with an
un-normalized shifter vector `d_i` (equivalently, shifting expenditure without
a counterparty budget).

## Why it seemed promising

A compositional shift is the natural way to express "households want more
housing"; it needs no new institutional block and leaves the production core
untouched.

## How it failed

The shifted household spends ≈ 25.5% more than its income
(`Σ_i β_i d_i p_i^{1−σ} ≈ 1.25485` at the reported prices). Household
expenditure exhaustion `Σ_i p_i c_i^h = E^h` fails, Walras' law can no longer
recover an omitted market equation, and the equilibrium identity chain breaks.
The failure was compounded by DE-0002.

## What replaced it

Two distinct, budget-complete experiments:

- **F1 preference reallocation** — renormalized weights
  `β̃_i = β_i d_i / Σ_j β_j d_j`, so exhaustion holds by construction
  (`cbase2/src/financing.jl`);
- **F2/F3** — a separate public-investment bundle `g_i` with an explicit
  counterparty `T(p)` or `F` (see ADR-0002, `ROADMAP.md` §4.3).

## Revival conditions

Only as F1 with `Σ_i p_i c_i^h = E^h` asserted exactly post-solve, or as an
explicitly financed injection. It must never be described as an autonomous
investment multiplier. Any revival also needs the household-expenditure
exhaustion gate from `ROADMAP.md` Phase 4.
