---
title: "IDEA-0009: pattern versus level in the sectoral response"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-10-24"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, sectoral-response, closure, ranking, referee-2, type-i-type-ii]
last-updated: October 2026
---

**Version 1** (October 2026)

Status: `scoped` (measured on the current runs; nothing in the manuscript, no ADR).

# The problem it addresses

The letter's answer to referee 2 claimed that "which sectors gain and which lose
is largely fixed by the input-output structure and the composition of the
program, and that ranking survives every closure". The manuscript makes no such
claim, and Section 5.2 shows only three sectors across the closures. The runs do
carry the sectoral solutions, so the claim can be tested rather than asserted
--- and the test says something sharper than either version.

# What the runs say

Data basis: `runs/matrix_5x3-v7|v8|v9-*` carry `solution.csv` for all fifteen
cells and match the reported `matrix_5x3-v11-*` runs exactly on all four
headline metrics (maximum deviation 0.0000; v8 and v9 bit-identical), so they
are the reported economies. Canary: $\sum_i (v_i/y_i)\,\Delta q_i$ equals the
manifest's $\Delta Y$ to $10^{-11}$.

- **The financing regime scales the response, it does not reorder it.** Within
  any labour closure, F1 versus F2 gives rank correlations of +0.97 (GAMMA,
  DELTA) to +1.00 (BF, ALPHA, BETA) across all 71 sectors, even where the
  aggregate response scales by a factor of six (GAMMA: $|\sum \Delta q|$ =
  0.0009 versus 0.0054 of GDP).
- **Within a price regime the labour closure does not reorder it either.**
  Across BF/ALPHA/BETA, and within GAMMA/DELTA, rank correlations are +0.92 and
  above with 94 to 100 per cent sign agreement.
- **It reorders only where the closure lets the response run.** GAMMA/DELTA in
  F3 flip the sign for 54 of 71 sectors relative to the fixed-price closures,
  and never in the other direction (zero sectors flip from + to $-$). Mean
  $|\Delta q|$ is 2 to 5 $\times 10^{-4}$ of GDP per sector, with only 2 to 7
  sectors moving more than 0.1 per cent of GDP, so most flips are near-zero
  crossings rather than reallocations.
- **The mechanism is a level shift, not a factor-intensity or relative-price
  reshuffle.** Rank correlation of $\Delta q$ with the wage share is +0.05 in
  both regimes and with the sectoral price change $-0.03$. In the fixed-price
  closures the programme crowds out output elsewhere ($\sum \Delta q$ is
  0.001 of GDP); where the wage is flexible the pool grows and the same sectors
  expand ($\sum \Delta q$ = 0.038, i.e. 2.86$\times G_0$ of gross output
  against 1.30$\times G_0$ in the Type I benchmark).

# Why it may matter for the paper

It is the paper's own thesis with a measurable edge: **the input-output structure
fixes the pattern of the sectoral response; the labour-market closure fixes its
level, and with it whether the response is a reshuffle of a fixed total or an
expansion in which almost every sector grows.** It also gives the honest scoping
for any sentence about which sectors gain and which lose --- true of the pattern,
not of the signs, and only where the response is material.

# What a specification would cost

- A rank-correlation table over both dimensions (five labour closures crossed
  with three financing regimes; all 71 sectors and the top 20 by size of
  response), plus the transition matrix for the one unstable pair.
- Data: the v7/v8/v9 solution files, already committed. Cost: an afternoon of
  scripting, no new runs. The only open choice is where "material" starts (0.1
  per cent of GDP? the top 20?).

# Open decisions for the operator

- State the pattern/level split as a sentence in Section 5.2, or leave Section
  5.2 and keep it in the letter? (Recommendation: one sentence in 5.2 --- R1.7
  asks for exactly this kind of sectoral stress-testing.)
- Publish the rank-correlation table, or the sentence only?
- Does the Type I/Type II labelling of the import content (letter R2.4,
  Appendix C.9) belong here or in its own note? (Recommendation: its own --- it
  is about the measure, not about the sectoral pattern.)

# Revision Log

- **Version 1** (October 2026) --- recorded the measurements, the mechanism, the
  cost of a specification and the scoping consequences.
