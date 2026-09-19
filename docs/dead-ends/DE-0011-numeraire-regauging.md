# DE-0011 — Re-gauging the numeraire to give eta_s a channel

- **Status:** dead end — rejected with a derivation and a measurement
- **Recorded:** 2026-09-18
- **Origin:** the "restore variety inside BF/ALPHA/BETA" candidate list
  (session prompt, 2026-09-18)
- **Related:** ADR-0019 (all-N with F), ADR-0014 (real-wage supply, DE-0004),
  BETA, `experiments/probes/numeraire_invariance.jl`,
  `experiments/probes/supply_identification.jl`
- **Evidence:** `experiments/probes/numeraire_invariance.jl` (demand-only gauge,
  2026-09-18), `experiments/probes/supply_identification.jl` (supply-shock
  pilot, 2026-09-18), `experiments/probes/analytic_headlines.jl` (closed forms,
  2026-09-18)

## What was tried

The proposal started from the observed BETA/ALPHA coincidence: with the CPI as numeraire the real wage did not move under the demand-only programme, so L^s = Lbar * ((w/P)/(w0/P0))^eta_s returned Lbar for every eta_s and the BETA row collapsed onto ALPHA. The suggestion was to adopt a wage numeraire (or another deflator) so that a demand shock would move the real wage and the supply elasticity would start to bite.

## Why it seemed promising

The gauge does fix the level of w when CPI(p) = 1, so it is tempting to conclude that changing the gauge frees the real wage. That reading mistakes the level for the ratio; it is the initial hypothesis this record rejects.

## How it failed

The real wage is not a numeraire object — it is a **technology** object. The
zero-profit block `p_i = c_i(p, w)` is homogeneous of degree 1 in `(p, w)` and
contains no demand term, so it determines only the ratio `p/w = p*(A)`. The
numeraire fixes the *level*: with `CPI(p) = 1` it gives `w = 1/CPI(p*)`; with
`w = 1` pinned it gives `p = p*`. Either way `w/P = 1/CPI(p*)`, which at the
calibrated technology (`p* = 1`, i.e. no productivity shock) equals 1 for
**any** demand vector. Re-gauging cannot create a channel that the zero-profit
block does not contain; it only rescales `p` and `w` by the same factor.

Re-measured on the ADR-0019 kernel (`experiments/probes/numeraire_invariance.jl`,
mobile ALPHA, programme scaled by k = 1/2/5/10): demand-only scaling leaves
`w = 1.0000000000`, `L = 1.0` and `max|dp| <= 4.1e-15` at every k. A programme ten
times the size leaves `w`, `p` and `L` at their baseline values, so there is
nothing for a different numeraire to re-gauge.

The pre-ADR-0019 mobile closed forms are retired along the same measurement:
with `F` entering household income `E`, mobile F2/F3 `consumption_rel` coincide at `-1.82262e-2` (k = 1; welfare index `0.981773801129`) and scale linearly in k, while income-side
`gdp_rel = 0` in every mobile cell (ADR-0018). At the BF row (`eta = 0`,
`p = w = 1`, `F = 0`) the closed form holds exactly
(`experiments/probes/analytic_headlines.jl`).

## What replaced it

Technology shocks move `w/P`. The +20% sector-1 pilot
(`experiments/probes/supply_identification.jl`) gives `w = 1.0049029745` and
separates `eta_s = 0.5 / 2 / 5` at `L = 1.0024484897 / 1.0098299882 /
1.0247564458`, with the implied elasticities recovered exactly. Alternatives
that could make demand move the real wage are sector-specific wages at
`eta = 0` and demand-sensitive pricing (markup or capacity) — never the gauge.

## Revival conditions

None via the numeraire. Only if the price block gains a demand term (markup or
capacity) would a demand-sensitive real wage appear — and then the channel
would live in the pricing rule, not in the gauge.
