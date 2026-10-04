---
title: "Revision narrative and commented outline"
project: "BFRep/(3)BeyondHulten"
date: 2026-09-20
version: 3
status: "working note - not a decision, no registry entry"
tags: [revision, narrative, outline, motivation-mindmap, beyondhulten]
---

**Version 3** \textcolor{revisionV2}{(September 2026)}
**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

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

## The section plan

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

## Key related references \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The table below repeats the first two columns of the
section plan and adds a third column naming, for every section of the
proposed paper, the key related references, embedded as links and verified by
web search (publisher or DOI pages preferred; the first occurrence carries
the link, later rows point back to it to keep the table readable). Several
entries also live in the repository bibliography; the list is agent-suggested
and remains to be pruned by the author.}

| Sect. | Working title | Key related references |
|---|---|---|
| 1 | Introduction | [Kuhn (1962)](https://press.uchicago.edu/ucp/books/book/chicago/S/bo13179781.html): inter-paradigm cleavages. [Kapeller (2013)](https://doi.org/10.1017/s1744137413000052): Model-Platonism and axiomatic variation. [Akerlof (1970)](https://www.jstor.org/stable/1879431): the narrative-legitimation analogy carried over from the old draft. |
| 2 | From transformation research to IO models | [Hornykewycz et al. (2025)](https://doi.org/10.1038/s44168-025-00229-2): the green-investment impulse. [Creutzig et al. (2018)](https://www.nature.com/articles/s41558-018-0121-1): demand-side transformation research. [Miller and Blair](https://www.cambridge.org/core/books/inputoutput-analysis/431D85A5E459AB078479852169EA77D7): input-output analysis. [Leontief (1936)](https://doi.org/10.2307/1927837) and [Leontief (1986)](https://ideas.repec.org/b/oxp/obooks/9780195035278.html): the IO core. |
| 3 | Two primitive questions | [Taylor (1990)](https://eclass.uoa.gr/modules/document/file.php/ECON249/4.%20%CE%A0%CE%9A%CE%9B%20%CE%BA%CE%B1%CE%B9%20%CE%A5%CF%80%CE%BF%CE%B4.%20%CE%93%CE%B5%CE%BD.%20%CE%99%CF%83%CE%BF%CF%81%CF%81%CE%BF%CF%80%CE%B9%CE%B1%CF%82/taylor%20STructuralist%20CGEs.pdf): structuralist closures. [Lavoie (2014)](https://www.e-elgar.com/shop/gbp/post-keynesian-economics-9781847204837.html) and [Blecker (2019)](https://www.e-elgar.com/shop/usd/heterodox-macroeconomics-9781784718893.html): post-Keynesian and heterodox macroeconomics. [Kim (2017)](https://arxiv.org/abs/1608.01365) and [Klump (2012)](https://doi.org/10.1111/j.1467-6419.2012.00730.x): the CES middle ground. [Weber (2023)](https://doi.org/10.4337/roke.2023.02.05), [Weber (2024)](https://doi.org/10.1093/icc/dtad080) and [Nikiforos (2024)](https://doi.org/10.1093/icc/dtae003): the Leontief price-model strand. [BF (2019)](https://doi.org/10.3982/ecta15202): the CES-GE pole. [Robinson (2006)](https://doi.org/10.1007/0-387-29748-0_11): multipliers and macro models across traditions. [Rocchi et al. (FIDELIO 3 manual)](https://publications.jrc.ec.europa.eu/repository/handle/JRC115308) and the [JRC FIDELIO model page](https://joint-research-centre.ec.europa.eu/projects-and-activities/trade-and-industrial-policy-analysis/industrial-policy/fidelio-model_en): the FIDELIO alternative. |
| 4 | One kernel, five closures, three financings | [BF (2019)](https://doi.org/10.3982/ecta15202), [Baqaee and Farhi (2021)](https://doi.org/10.1257/pandp.20211107) and [BF (2022)](https://www.nber.org/papers/w27152): the closure lineage and the slack argument. [Destatis input-output accounts](https://www.destatis.de/DE/Themen/Wirtschaft/Volkswirtschaftliche-Gesamtrechnungen-Inlandsprodukt/Tabellen/_tabellen-innen-in-output.html): the A-bill calibration data (rows 73-75). |
| 5 | Application: the programme through the matrix | Hornykewycz et al. (2025; see the row for Section 2): impulse provenance and price bases. The sectoral-labour evidence base: [IAB-Stellenerhebung Q4/2024](https://iab.de/presseinfo/iab-stellenerhebung-fuer-das-vierte-quartal-2024-zahl-der-offenen-stellen-steigt-saisonbedingt-auf-14-millionen/), [KfW-ifo skills barometer, December 2024](https://www.kfw.de/PDF/Download-Center/Konzernthemen/Research/PDF-Dokumente-KfW-ifo-Fachkr%C3%A4ftebarometer/KfW-ifo-Fachkraeftebarometer_2024-12.pdf) and the [ifo order-book range report](https://allgemeinebauzeitung.de/abz/auftragsreichweite-im-bauhauptgewerbe-auftragsbestaende-ruecklaeufig-57948); details in `construction_capacity.tex`. |
| 6 | Beyond demand-only shocks | [Samuelson (1951)](https://cowles.yale.edu/research/cfm-13-activity-analysis-production-and-allocation): the nonsubstitution theorem behind Lemma 2. [Kaldor (1961)](https://doi.org/10.1007/978-1-349-08452-4_10) and [Sylos-Labini (1995)](https://doi.org/10.1016/0954-349x%2895%2900025-i): the demand-led productivity and production-function doors. BF (2019; see the row for Section 3): the frictionless reallocation pole. |
| 7 | Discussion | BF (2019, 2021, 2022; see the row for Section 4): production networks and the COVID-19 benchmark; plus [Acemoglu, Akdoglu and Kerr (2015)](https://www.nber.org/papers/w21344) on networks and the macroeconomy. Taylor (1990), Rocchi et al. and Weber (2023, 2024) / Nikiforos (2024) as linked under Section 3. |
| 8 | Conclusion | Kapeller (2013) and Kuhn (1962), linked under Section 1: the methodological frame. |

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

# External review \textcolor{revisionV2}{\normalsize [added v3]}

\textcolor{revisionV2}{This chapter records the agent's read-only currency
check of the outline against repository HEAD \texttt{e78ae88} (2026-09-20,
working tree clean), carried out on 2026-09-30 immediately before the
drafting pass was to begin. It is a report, not a decision: the body of this
note above is unchanged from v2, and none of the findings below has been
applied. The author decides which to take up.}

## \textcolor{revisionV2}{Verdict}

\textcolor{revisionV2}{The outline is current as a vantage point but not yet
draft-ready. Nothing in the repository has moved since it was written --- the
outline itself was last touched in \texttt{6026806} (improved narrative) and
no later commit touches it --- and every source it cites exists at the stated
version. What does not survive contact with the repository is (i) three of
the outline's own factual claims, one of which is a substantive equivalence
error, and (ii) two older manuscript artefacts it silently collides with: the
front matter of \texttt{revised\_manuscript/} and the Stage 3 plan in
\texttt{docs/DOCS\_ASSESSMENT.md}, now at Version 8.}

## \textcolor{revisionV2}{What checks out}

\textcolor{revisionV2}{The table below is new in v3; it is left uncoloured
under the no-coloured-markdown-table rule.}

| Claim in the outline | State at HEAD |
|---|---|
| `equivalence.tex` v4, evidence layer `matrix_5x3_v10` (C1-corrected, 33/33) | Confirmed in the file header (commit `ecdf814`); v10 is the citable generation for the matrix and the eighteen sectoral cells, v9 for everything else |
| Supply arm executed, ETAs v5 | `ETAs.md` is v5 (revisionV4); **189** cells, not 126: `supply_etas_s1` (63, the identification vehicle), `supply_etas_prog` (63) and `supply_etas_unif` (63) |
| Sectoral family ADR-0022 on `matrix_5x3_v9/v10` | Confirmed; v9 superseded only for the eighteen sectoral gdp-family metrics (defect C1), and the v9 runs stand as issued (ADR-0004) |
| Reviewer maps R1.1-R1.9, R2.1-R2.14 | Present in `docs/DOCS_ASSESSMENT.md`, which is at **Version 8** |
| Ladder monotone in rigidity; F2 max abs(p-1) 0.153 to 0.037 | Confirmed on v10: `BETA-F2-etas{025,05,1,2}` = 0.152897 / 0.105666 / 0.065422 / 0.037171 |
| Grouping rule not yet ratified | Confirmed: `docs/grouping_rule_evidence.md` v1 still states the rule is for ratification |

## \textcolor{revisionV2}{Corrections required}

### \textcolor{revisionV2}{The F2 equivalent-of-F3 shading in the fixed-wage row}

\textcolor{revisionV2}{The presentation rule crosses the three labour classes
with (F1, F2 $\equiv$ F3). That is wrong in the GAMMA/DELTA row: financing
neutrality holds only where the external transfer $F$ is the endogenous
instrument --- the BF row at $\eta = 0$, ALPHA/BETA, and all six sectoral
variants, where F2 and F3 agree to all printed decimals --- and fails exactly
in the fixed-wage rows, whose system pins $F \equiv 0$. Measured on the v10
manifests:}

| Cell | Consumption rel. | Employment | External position |
|---|---:|---:|---:|
| `GAMMA-F2` | -1.533 % | 1.00126 | 0.0 |
| `GAMMA-F3` | **+2.264 %** | **1.01780** | +0.01331 |

\textcolor{revisionV2}{\texttt{DELTA-F2} and \texttt{DELTA-F3} show the same
split, and \texttt{paper/equivalence.tex} already says so (the F-instrument
paragraph and the "who pays" paragraph: F3 expands employment by 1.78 \%
and consumption by 2.26 \%, the Keynesian case, against a tax-financed
welfare loss under F2). This is not only a shading fix: the fixed-wage row is
where the financing choice actually bites, so it belongs in Section 5's
robust-versus-sensitive table. Consequence for the abstract claim: the
fifteen cells collapse to seven distinct economies --- BF two, ALPHA
$\equiv$ BETA two, GAMMA $\equiv$ DELTA three --- not to three; three is a
count of labour classes.}

### \textcolor{revisionV2}{The GDP share of the programme}

\textcolor{revisionV2}{The narrative arc states that G0 = 40.3 bn in 2019
prices equals the horizon mean, 1.945 percent of GDP. That sentence mixes
price bases: 1.945 \% is the raw 2023-price mean over the 2019-basis
denominator (\texttt{GDP\_P} = 3,027,818 EUR m; 58.893 / 3027.818). The
model's programme share is 1.331 \% (40.3 / 3027.8, equal to the booked
\texttt{B\_gov} = 0.01330991), and ADR-0024 attaches the about-1.9 \%-of-GDP
figure to the 2024 first-year peak (81 bn, 2023 prices), not to the mean.}

### \textcolor{revisionV2}{The G0 open item}

\textcolor{revisionV2}{The open item names the ADR-0023 fork on G0 and asks
for it to be settled before Section 5's tables. The fork is ADR-0024, and it
is accepted (2026-09-20): keep G0 = 40,300, narrative-only, no re-anchor, no
re-run, no re-mint. ADR-0023 is the supply-shock design schema (accepted,
implemented, executed). What remains open is the composition of psi --- the
raw 2024 current-price incidence vector; a constant-2019-price vector needs
sectoral deflators and would move the incidence and rigid-group results ---
and the time-aggregation note.}

### \textcolor{revisionV2}{Pointer and coverage defects}

- \textcolor{revisionV2}{The source list undercounts the supply arm: 189
  cells across three designs, not 126 (the table above).}
- \textcolor{revisionV2}{Appendix C should name the tally it holds: fifteen
  matrix cells plus twelve ladder cells plus six two-group cells = 33.}
- \textcolor{revisionV2}{\texttt{paper/framing\_gamma.tex}, the
  manuscript-facing GAMMA note, is referenced nowhere in this outline; and
  the ADR-0021 wage-structure axis (still proposed) is a candidate
  representation-robustness row the presentation rule does not mention.}

## \textcolor{revisionV2}{Factors to consider when setting out}

- \textcolor{revisionV2}{Where the draft lives --- decide first.
  \texttt{revised\_manuscript/} already exists: \texttt{main.tex} (29 July)
  plus eight chapter files, all marked as old text to be revised under a
  salvage plan, with the new per-cell panels in \texttt{pictures/}. Its front
  matter carries the retired thesis: the labour-supply elasticity as the key
  parameter, and the claim that the factor-market closure explains roughly an
  order of magnitude more variation than all production technology parameters
  combined. The second sentence is on the Stage 3 must-not-appear list. Either
  the new draft replaces that front matter in place, keeping ch01 and ch02 as
  the salvage, or it lives in a fresh file; Sections 3 to 6 of the outline
  have no chapter homes yet.}
- \textcolor{revisionV2}{Two competing tables of contents. The Stage 3 plan
  in \texttt{docs/DOCS\_ASSESSMENT.md} (v3-era, revisionV2 colour) proposes a
  nine-section manuscript built on the two wedges, with sections Aggregate
  results, Sectoral results and The financing wedge, and it holds the reviewer
  point maps --- the one place where every referee item has a named location,
  from which the response letter is built. This outline re-homes the narrative
  but leaves those maps unanchored. Recommend a reviewer-coverage column in
  the section plan, and importing the Stage 3 must-not-appear list (all
  eta-dominance claims, the 88.4 \% and 100 \% figures, GO certification
  language, price-invariance claims, and the claim that the IO multiplier
  requires eta to infinity) into the retirement section.}
- \textcolor{revisionV2}{Build path. The draft cannot be compiled in the
  container: \texttt{latexmk -lualatex} stops at a missing
  \texttt{authblk.sty}, and the Nix texlive environment, 542 packages, has no
  biblatex, biber, pgfplots or makecell either. The manuscript is a Mac
  compile artefact today; decide before building a write-compile-verify loop
  on it --- keep Mac compilation, or extend the container toolchain, which is
  a large download and needs the author's call.}
- \textcolor{revisionV2}{The repository contract is stale.
  \texttt{AGENTS.md}, last commit \texttt{ac49c5d} of 19 September, still
  names \texttt{matrix\_5x3\_v5} as the generation to cite, reports
  \texttt{DOCS\_ASSESSMENT.md} as Version 5, and knows nothing of ADR-0020 to
  ADR-0025, the sectoral family or the supply arm. Refresh it before drafting
  so the invariants list matches the generations being cited; it is a
  protected file, so the write needs an interactive approval.}
- \textcolor{revisionV2}{Unratified dependence. Section 6's rigid-group cells
  rest on the grouping rule, which is not ratified; today only the two
  executed rules are citable (programme sectors, 0.272033; largest half,
  0.123362). Either ratify first or write Section 6 to the executed rules.}
- \textcolor{revisionV2}{Open Stage 0 items that bite the draft rather than
  the model: the audit-trail disclosure level, and the scope locks (skill
  classes, Armington, the B\&F replication off the critical path).}
- \textcolor{revisionV2}{Measurement semantics: real GDP is the income-side
  index, flat by construction under demand-only shocks (ADR-0018), while the
  consumption index is welfare. The old abstract's claim that the GDP effect
  ranges from virtually zero to about 1.5 \% conflates the two, and that
  sentence belongs to the retired framing.}
- \textcolor{revisionV2}{Evidence discipline worth advertising: every table
  cites run ids, the v9/v10 citation split matters, and the preregistration
  trail with 41 failed and 3 provisional runs left visible is a credibility
  asset for this submission. The table of contents has no home for it; a
  short reproducibility subsection in Appendix B would carry it.}

## \textcolor{revisionV2}{Suggested edits}

1. \textcolor{revisionV2}{Rewrite the presentation rule: columns F1, F2, F3;
   shade F2 $\equiv$ F3 in the BF and ALPHA/BETA rows only; state the
   neutrality condition, that $F$ is endogenous; flag the fixed-wage row as
   the financing-sensitive case; restate the collapse as fifteen cells to
   seven.}
2. \textcolor{revisionV2}{Fix the Section 1 numbers: 40.3 bn per year at 2019
   prices, about 1.331 \% of GDP; 58.9 bn per year at 2023 prices; 2024 peak
   81 bn, about 1.9 \% of GDP.}
3. \textcolor{revisionV2}{Close the G0 open item (ADR-0024 accepted, option
   A) and replace it with the composition and time-aggregation items; correct
   the ADR reference.}
4. \textcolor{revisionV2}{Restate the supply arm as 189 cells, naming the
   three designs.}
5. \textcolor{revisionV2}{Add the 33-cell tally and the ADR-0021
   wage-structure axis to Appendix C, and reference
   \texttt{paper/framing\_gamma.tex} in Sections 4 and 6.}
6. \textcolor{revisionV2}{Add the reviewer-coverage map and the
   must-not-appear list.}
7. \textcolor{revisionV2}{Add the reproducibility item to Appendix B, and
   the draft-location and build-path decisions as named open items.}

## \textcolor{revisionV2}{Decisions still owed}

- \textcolor{revisionV2}{Is the new draft written into
  \texttt{revised\_manuscript/chapters/}, reusing ch01 and ch02 and retiring
  the rest, or as a fresh main file?}
- \textcolor{revisionV2}{Is \texttt{AGENTS.md} refreshed in the same pass? It
  is protected, so it needs an interactive approval.}

# Revision Log

- **Version 1** (September 2026)
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Section 5
  restructured: the section plan is now 5.1, and a new 5.2 reproduces its
  first two columns and adds a third column with the key related references
  per section, web-verified and embedded as links (same-round addition,
  folded into v2).
- **Version 3** \textcolor{revisionV2}{(September 2026)} --- new final
  chapter "External review" (placed before this log so the log stays last):
  the agent's read-only currency check of the outline against repository HEAD
  `e78ae88` (2026-09-20), with the pointer claims that check out, the
  corrections required (the F2 $\equiv$ F3 shading in the fixed-wage row, the
  GDP-share attribution, the settled G0 item), the factors to weigh before
  drafting, the suggested edits and the two decisions still owed. The body
  above is unchanged from v2 and no suggested edit has been applied.