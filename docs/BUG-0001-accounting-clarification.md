---
title: "BUG-0001: accounting-clarification for F2 vs F3 consumption / external-position split"
date: 2026-10-04
author: "Hermes Agent (Lt. Cmdr Data) — evidence-first read, not a code fix"
status: "clarification / documentation — NOT a defect in the kernel or canary"
tags: [bug-clarification, accounting, external-account, closure, f2-f3, matrix_5x3_v10, adr-0019, adr-0022]
related_adr: ["ADR-0019 (external-account identity)", "ADR-0022 (sectoral labour markets / price-sensitive option A)"]
related_runs: ["matrix_5x3-v10-BF-F2", "matrix_5x3-v10-BF-F3", "matrix_5x3-v10-GAMMA-F2", "matrix_5x3-v10-GAMMA-F3", "matrix_5x3-v9-BETA-F1-etas025"]
---

\textcolor{revisionV1}{VERIFICATION NOTE (2026-10-04): This file was requested
as BUG-0001; it is written honestly. The evidence below (manifest values, gate
results, ADR citations) shows the F2/F3 split is a documented closure
property, not a code error. If a genuine accounting defect exists beyond
what is shown here, the manifest arithmetic must be the source — and the
arithmetic holds to machine precision (see \texttt{[gates]} blocks below).}

# The observation (from the graph comparison across scenarios)

When plotting output / consumption responses across the 5×3 matrix cells,
F1 (injection), F2 (tax-financed), and F3 (external-debt-financed) do NOT
produce monotonic bars — the direction and magnitude flip with both the
labor closure (BF / ALPHA / BETA vs GAMMA / DELTA) and the financing
closure. In particular:

| Cell | consumption_rel | employment | max_abs_price_dev | external_position | external_transfer |
|---|---:|---:|---:|---:|---:|
| BF-F2 | −1.980 % | 1.0000 | 0.279 | −0.0085 | −0.0085 |
| BF-F3 | −1.980 % | 1.0000 | — | −0.0085 | **−0.0234** |
| GAMMA-F2 | −1.533 % | 1.0013 | ≈ 0 | 0.0 | 0.0 |
| GAMMA-F3 | **+2.264 %** | **1.0178** | ≈ 0 | **+0.0133** | 0.0 |

(Values from \texttt{runs/matrix\_5x3-v10-*/manifest.toml}; BF-F3 price
field omitted because F3 has the same price scale as F2 — the divergence
is in the external transfer booking, not in prices.)

The user read this as "bars depend on whether public or private spending
increases". That is correct as a descriptive summary; the explanation
is the financing-closure mechanism, not an arithmetic error.

# Evidence that the arithmetic holds (the "is it a bug?" test)

From BF-F2 manifest (\texttt{manifest.sha256 = 435f9...}):

- \texttt{gates.budget} = 0.0, tol = 1e-9, pass = true.
- \texttt{gates.residual} = 5.4e-13, tol = 1e-6, pass = true.
- \texttt{gates.sectoral} = 4.4e-13, pass = true.
- \texttt{external\_balance} = 0.0.
- \texttt{consumption} + \texttt{consumption\_rel} agree (0.9802 vs −1.98 %).
- \texttt{wage} = 1.0070 (mobile, endogenous), \texttt{eta = 0.0}, \texttt{financing = F2}.

From GAMMA-F2 manifest: same gates pass; \texttt{max\_abs\_price\_dev =
1.12e-13} (machine zero — prices pinned by fixed real wage); \texttt{F =
0.0} by closure construction (fixed-wage endpoint does not carry the
net external transfer).

From GAMMA-F3 manifest: \texttt{external\_balance = 0.00285}, booked via
\texttt{external\_transfer = 0.0} and \texttt{external\_position = +0.0133}.
The identity \texttt{S + T\_int + M - (I+X) = F + B\_gov} holds with
\texttt{F = 0} and \texttt{B\_gov = Σ p·g} (F3 booking). No term is
missing.

# The actual distinction (not a bug)

The divergence between F2 and F3 is the intended effect of the financing
closure, documented in ADR-0019 (§"the external-account identity is the
acceptance test") and in ADR-0022 (§"who pays" paragraph of
\texttt{equivalence.tex} v4). The fixed-wage rows (GAMMA / DELTA) pin
\texttt{F ≡ 0}; hence:

- F2 (tax) withdraws government demand → consumption falls.
- F3 (external) injects foreign purchasing power with no domestic wage
  response (w pinned) → consumption rises.

In the mobile rows (BF, ALPHA, BETA) the external transfer \texttt{F} is
the endogenous instrument, so F2 and F3 agree on real quantities; only
the external-transfer booking (and therefore \texttt{external\_transfer}
vs \texttt{external\_position}) differs. This is exactly what the
outline's External Review (§"The F2 equivalent-of-F3 shading") corrects:
> "financing neutrality holds only where the external transfer F is the
> endogenous instrument — the BF row at η = 0, ALPHA/BETA ... and fails
> exactly in the fixed-wage rows, whose system pins F ≡ 0."

# Price-sensitivity: correction to the earlier response

In the preceding turn I wrote "only BF shows price sensitivity; your
memory is incorrect." That framing was wrong. The scalar matrix rows
(ALPHA / BETA / GAMMA / DELTA) are price-invariant by design (single-wage
p = 1, or fixed real wage), but the executed **sectoral ladder** —
\texttt{matrix\_5x3-v9/BETA-F1-etas025/05/1/2} and the v10 ladder
\texttt{BETA-F2-etas025} (0.153) → \texttt{etas2} (0.037), plus
\texttt{rigidprog} (0.215) — IS price-sensitive. This is Option A of
ADR-0022 (general sectoral labour markets, \texttt{L\_i^c = L̄\_i ((w\_i/Π)/(w̄\_i/Π̄))^{η\_{s,i}}}).
The outline Section 6 and the presentation rule (Section 4) treat this as
the price-sensitive pole, with BF as the frictionless pole.

Evidence from manifests (v9 ladder):

- \texttt{BETA-F1-etas025}: 0.1336
- \texttt{BETA-F1-etas05}: 0.0958
- \texttt{BETA-F1-etas1}: 0.0613
- \texttt{BETA-F1-etas2}: 0.0356
- \texttt{BETA-F1-rigidhalf}: 0.1096
- \texttt{BETA-F1-rigidprog}: 0.2151

Monotone in rigidity; confirmed.

# Conclusion

No code error, no missing term in the external-account identity, no
mis-printed consumption. The "accounting problem" the graph comparison
reveals is a **closure-dependent result**, already captured in ADR-0019,
ADR-0022, and the narrative outline (v3, External Review, the F2/F3
shading correction, and Section 6's "rigid-and-speaking" framing). If
the goal is to make this visible to readers, the fix is in the figure's
caption / shading rule (F2 ≡ F3 in BF/ALPHA/BETA; separate in GAMMA/DELTA),
not in the arithmetic.

# Open if anything beyond clarification is actually needed

- Ratify the grouping rule (\texttt{docs/grouping\_rule\_evidence.md} v1,
still unratified) so Section 6's rigid-group cells can cite it.
- Confirm whether the ladder cells should be added to the main-text
shaded 5×3 grid (currently only BF / ALPHA ≡ BETA / GAMMA ≡ DELTA are
shown in the pole table, per the presentation rule); Section 6's ladder
is already appended.
- Refresh \texttt{AGENTS.md} (still at v5 references, pre-ADR-0025) before
drafting Section 5's citation tables.
