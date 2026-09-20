---
title: "Revision narrative and commented outline"
project: "BFRep/(3)BeyondHulten"
date: 2026-09-20
version: 1
status: "working note - not a decision, no registry entry"
tags: [revision, narrative, outline, motivation-mindmap, beyondhulten]
---

# Source and status

This note consolidates the narrative and table-of-contents proposal for the
revised paper. It is derived from:

- `docs/motivation_mindmap.md` (v3) - the author-confirmed semantics of the
  MOTIVATION board;
- `paper/equivalence.tex` (v4, evidence layer `matrix_5x3_v10`, the
  C1-corrected 33/33 generation) - the ex-ante equivalences and their scope;
- the executed supply-side arm (`supply_etas_prog|unif`, 126 cells; ETAs v5)
  and the sectoral-labour family (ADR-0022, `matrix_5x3_v9/v10`);
- the reviewer maps in `docs/DOCS_ASSESSMENT.md` (R1.1-R1.9, R2.1-R2.14);
- four framing decisions put to and taken by the author on 2026-09-20
  (next section).

Everything not covered by those decisions is agent-suggested wording:
propose-first, no decision weight. The note lives in `paper/` as a supporting
manuscript note (cf. `closures_and_demand_shocks.md`, `framing_gamma.tex`);
it is not a plan file and carries no registry entry.

# Decisions taken (2026-09-20)

| Open node | Decision |
|---|---|
| Headline result | The transformation use case first; the equivalences serve as supporting rigour. Do not say "collapse theorem" in the abstract. |
| `BF?` | Frictionless pole and boundary case, implemented with peer-modesty: introduced late in the framing, no vantage status. |
| PARADIGMS hunch | Confirmed: B and C foreground; A derived conceptually but becomes the measured object from the matrix chapter on; D and E in passing. |
| `WHICH GAMMA? / WHICH BETA?` | A 3x3 pole table in the main text with equivalence shading; the full 15-cell matrix and all gamma/beta variants in the appendix. |

# Thesis and working title

**Working title:** *Which Short Run? Which Money? Closure Choice in
Input-Output Models of the Socio-Ecological Transformation.*

**Thesis sentence:** the transformation short run is doubly rigid ---
quantities do not arbitrage (the IO heritage) and sectoral wages do not
clear (the sectoral-labour family) --- and BF is the frictionless pole
against which that rigidity is measured.

**Spine in one sentence:** applied transformation research needs sectoral IO
models, but every IO calculation silently embeds two macro commitments --- a
vision of the short run (B) and a theory of money (C) --- whose joint
operationalisation is the adjustment closure (A); we map that space as a 5x3
closure matrix run on an empirically informed green-investment vector, show
that under demand-only shocks the matrix reduces to a few distinct economies,
and identify exactly where and why it does not.

The one refinement worth defending against the plain "more rigid" reading:
relative to BF our closures are indeed more rigid, and the executed ladder is
monotone in rigidity (F2 max abs(p-1) falls from 0.153 to 0.037 as
$\eta_s$ rises from 0.25 to 2); but relative to the pure Leontief endpoint
they are *less* rigid on the price margin --- DELTA under demand-only shocks
is rigid-and-silent ($p \equiv 1$ exactly), the sectoral family is
rigid-and-speaking. The contribution is the map of rigidities, not a
monotone rigidity claim.

# Narrative arc

1. **Concrete, then abstract.** Open with the transformation programme (the
   npj Climate Action impulse: about 58 bn/yr and 3.1 tn until 2050 in 2023
   prices; G0 = 40.3 bn in 2019 prices = the horizon mean, 1.945 percent of
   GDP; cite both price bases per ADR-0024). Beginners grasp this
   immediately; the trap comes next.
2. **Name the hidden choice.** The "multiplier" everyone quotes assumes a
   short run and a financing regime. Closure is *the* axiomatic variation of
   this literature; the Kuhnian framing of the old Section 1 survives the
   revision intact.
3. **Operationalise.** The 5x3 matrix: five labour closures
   (BF, ALPHA, BETA, GAMMA, DELTA) crossed with three financing closures
   (F1 preference reallocation, F2 tax-financed, F3 external debt). D and E
   appear here in passing: one fixed kernel for all cells, prices
   demand-invariant under demand-only shocks.
4. **Application.** Run the programme through the matrix. Which conclusions
   are robust, which are closure-sensitive --- this is the headline the
   reader is promised in the introduction.
5. **Rigour engine.** The equivalences (Lemma 1 price-block
   demand-invariance; Lemma 2 nonsubstitution; Propositions 1-3; financing
   neutrality) justify showing three economies instead of fifteen cells, and
   state precisely where closure choice is redundant. The no-go corollary
   proves the rigidity that matters cannot be manufactured anywhere else
   (nominal wage rules are vacuous in both gauges): it must live in
   sectoral labour markets.
6. **Relevance return.** Re-read the transformation application through the
   map: no closure-free multiplier exists, and the paper can now say *which*
   choice matters *for what* question.

# Commented table of contents

| Sect. | Working title | What the argument does | Mindmap |
|---|---|---|---|
| 1 | Introduction | Transformation hook; the hidden closure choice; contributions (rigidity map, closure taxonomy, 5x3 on a published impulse, equivalences listed but not headlined); reading-guide box | N01, N08-N10 |
| 2 | From transformation research to IO models | Why *plural* IO models are the bridge; construction and price basis of the empirically informed vector; light mathematics, the beginner on-ramp | N08, N09, N10 |
| 3 | Two primitive questions | 3.1 vision of the short run: quantity adjustment vs. rationing vs. clearing, mapped via Taylor closures, FIDELIO and the neo-Keynesian rigidity catalogue. 3.2 theory of money: F1/F2/F3 as endogenous-money vs. balanced-budget vs. external-funds answers (F3 claimed as an *accounting* closure, not a monetary mechanism). A stated as the operational summary of 3.1 and 3.2 | N02B, N02C, N02A |
| 4 | One kernel, five closures, three financings | The common kernel (E in passing: one production core deliberately held fixed); row and cell definitions; the "why only three economies" box (Propositions 1-3 in one paragraph, forward reference to Section 6 and Appendix A); the presentation rule as one paragraph | N02A, N06, N11, N12 |
| 5 | Application: the programme through the matrix | Multipliers, price responses, external position across cells; the robust-vs-closure-sensitive table; decision heuristics for transformation modellers | N06, N05, N10 |
| 6 | Beyond demand-only shocks | What simplifies and where variety lives: the equivalences and their scope (demand-only, the A = 1 branch the numeraire selects); the sectoral family and its exact nesting (the rigid corner *is* the $\eta = 0$ endpoint); the executed ladder and rigid-group cells; the supply-shock arm as identification vehicle; the no-go corollary; BF repositioned here as the frictionless pole | N03, N13, N02D |
| 7 | Discussion | Positioning against BF and the production-network literature (pole, not foil); FIDELIO and structuralist macro; limitations; open doors (capacity, markup) acknowledged as future work | N03, N02D |
| 8 | Conclusion | Return to the thesis: closure choice is first-order, and the matrix tells you where it is and is not | --- |
| App. A | Full proofs | Lemmas 1-2, Propositions 1-3, from `equivalence.tex` v4 | N13 |
| App. B | Kernel, calibration, A-bill accounting | Destatis rows 73-75, the external-account canary identity | N05 |
| App. C | Matrix variants and gamma/beta representation robustness | Where the presentation-decision variants live; all tables cite run ids | N11, N12 |
| App. D | Vector construction and data provenance | `impulses.csv`, the ADR-0024 price-basis verification | N08 |

Audience devices: the reading-guide box after the introduction (beginners:
Sections 1, 2 and 5; experts: Sections 3, 4 and 6); one figure doing the
heavy lifting --- the shaded 5x3 grid (candidates catalogued in
`paper/possible_plots.md`, ADR-0025); all proofs in Appendix A; all result
tables cite `run_id`s from the `matrix_5x3_v10` and `supply_etas_*`
generations.

# The 5x3 presentation rule

The equivalences dissolve most of the gamma/beta version question. Under
demand-only shocks ALPHA $\equiv$ BETA and GAMMA $\equiv$ DELTA are exact, so
"which version to show" in those rows is a choice of labels, not numbers; the
labels separate only where the equivalence fails --- the sectoral ladder and
the supply-shock arm --- and those cells are shown anyway in Section 6.

Concretely:

- **Main text, Table 1:** a 3x3 pole summary --- BF's sectoral-wage system;
  the single-wage flexible class (ALPHA, BETA); the fixed-real-wage quantity
  class (GAMMA, DELTA) --- crossed with (F1, F2 $\equiv$ F3), with the
  equivalence classes shaded.
- **Appendix C:** the full 15-cell matrix (citing `matrix_5x3_v10` run ids)
  plus the gamma/beta variant tables.
- **One figure:** the shaded 5x3 grid.
- Which beta to foreground: beta's substantive content in this paper is the
  sectoral family (ADR-0022) --- its canonical presentation is the executed
  ladder and the rigid-group cells, not the scalar-elasticity variant that
  collapses onto ALPHA.

# What the revision retires

Stating this explicitly in the manuscript is the honest structural answer to
the reviewer criticism:

- The exogenous endowment shift / unfinanced-slack mechanism of the old
  Section 5 (the +19.3 pp headline) is dropped and inadmissible: every shock
  is now financed (F1/F2/F3).
- The old elasticity-gradient and slack-calibration sections (Sections 4-5 of
  `paper/main.tex`) are superseded by the 5x3 matrix.
- External-position claims from superseded generations (the pre-ADR-0020 BF
  pin) are never printed; from `matrix_5x3_v6` on the BF economy is
  identified.
- Surviving from the old draft, largely as-is: Section 1 (Kuhnian cleavages,
  axiomatic variation) and Section 2 (IO-table representation).

# Pitfalls and guards

- Do not let A-E become five coequal sections: the 5x3 matrix *is* the
  decomposition of A into B and C; the taxonomy and the nesting are the same
  object. D is inert in the demand-only scope ($p \equiv 1$ in the
  single-wage rows) and E is held fixed by design --- in passing is not a
  preference but a consequence.
- Never claim bottleneck inflation from GAMMA under demand-only shocks: the
  measured price response is exactly 0.000. Price movement enters only via
  BF's sectoral wages or the sectoral family.
- The elastic end of the family is a limit, not a cell: the solve fails at
  $\eta_s = 10^6$, so a family cell coincides with GAMMA only asymptotically.
- Scope the equivalences airtight (demand-only, the A = 1 branch the
  numeraire selects); give a clean counterexample for supply shocks and use
  it to motivate Section 6 rather than hide the boundary.
- Keep the gamma/beta cluster as the one-paragraph representation rule in
  Section 4; leaking it as a debate reads as indecision.
- Fixed-wage rows have a documented warm-start sensitivity (about
  $1.6 \cdot 10^{-11}$): compare generations at $10^{-10}$, never by
  bit-identity.

# Open items

- The BF peer-modesty implementation (late introduction, no vantage status)
  is the agent's suggestion for executing the author's A-or-C answer; it
  needs the author's reading pass.
- The reading-guide wording and the choice of the one heavy-lifting figure
  among the `possible_plots.md` candidates.
- The ADR-0023 fork on G0 (narrative-only with a linear-scaling note vs.
  re-anchoring) interacts with Section 2's price-basis paragraph and should
  be settled before drafting Section 5's tables.
- Section 6's ladder presentation should be checked against the ratification
  status of the grouping rule (`docs/grouping_rule_evidence.md`), which
  determines the rigid-group cells to be cited.
