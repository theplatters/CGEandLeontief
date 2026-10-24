---
title: "Shortening candidates: the nineteen footnotes"
author: Hermes Agent
date: 2026-10-24
last-updated: October 2026
tags: [shortening, footnotes, cleaner-pass, beyondhulten, revision]
project: BeyondHulten
---

**Version 1** (October 2026)

This note collects the places where the manuscript can be shortened without
touching an argument, a number or a claim. It implements rule 8 of
`docs/prompts/manuscript-cleaner.md` ("Prefer none; where a footnote merely
restates a clause, fold it into the sentence or delete it") and is the
authorised home for the items the author has ruled out, so that they are not
re-opened on the next pass.

The inventory below was taken from the committed tree at `cf53520`, live text
only (LaTeX comments excluded), and counts words in the footnote body.

# Summary

The manuscript carries **19 footnotes, 1,068 words** --- roughly two pages of
small type across 68 pages. Ten of them are candidates; together they hold
**659 words**, 62 % of the footnote mass. The remaining nine are definitional
and should stay.

| Verdict | Footnotes | Words | Share |
| --- | --- | --- | --- |
| Relocate | 4 | 328 | 31 % |
| Fold into the sentence | 3 | 108 | 10 % |
| Cut | 3 | 223 | 21 % |
| Keep | 9 | 409 | 38 % |

Nothing here is a correctness fix. Every candidate is a formatting decision, so
the manuscript compiles and reads correctly if none of them is adopted.

# Relocate

A footnote that carries a paragraph of content is a paragraph in the wrong
place. All four below either support a claim the text already makes or state a
mechanism the appendices already derive.

- **§2 L9 (46 words) --- the sectoral incidence of the program.** Lists the
  seven shocked sectors with their expected relative changes (+50.3 % glass and
  glassware down to +1.0 % electrical equipment). This is the only place in the
  manuscript where the incidence is stated at that level of detail: the
  composition appendix (C.9) decomposes the bill into domestic, imported and
  taxed components but not by sector, and Section 5 repeats only three of the
  seven shares in prose. **Action:** add the seven-sector incidence as a panel
  of the composition table in Appendix C and reduce the footnote to a pointer.
- **§3 L75 (129 words) --- how FIDELIO and E3ME are wired.** The wage curve,
  mark-up pricing with sluggish pass-through, and the error-correction
  consumption rule. The text needs one claim from it (both models pin the two
  sides of the matrix); the mechanics are a second contribution to the
  comparison and belong beside the equivalences. **Action:** move to
  Appendix A and keep one sentence in §3.
- **§5 L28 (70 words) --- why $W_{CES}$ and $W_C$ differ.** The $\sigma=1$
  Cobb--Douglas equivalence. This is measure definition, and the measure is
  already discussed in Section 4.3 and tabulated in Table 5. **Action:** move
  beside the measure definition or into the table note.
- **§6 L15 (83 words) --- the continuation protocol.** Four starting points,
  agreement to $10^{-6}$, full-rank residual Jacobian. This is replication
  detail, not an economic qualification. **Action:** move to the replication
  package or the uniqueness guard in Appendix E.

# Fold into the sentence

Three footnotes restate a clause that is already on the page. Each can be
absorbed in a few words, which is what rule 8 asks for.

- **§3 L47 (48 words).** "Assumptions deviate from the neoclassical baseline but
  strengthen the general laissez-faire policy stance" plus a signpost citation.
  **Action:** append the citation to the sentence and drop the rest.
- **(A) L31 (28 words).** "Considerations from this appendix do not necessarily
  apply there." **Action:** fold into the preceding sentence of Appendix A.
- **(B) L537 (32 words).** Writing $\eta_{s,i}$ for $\eta_s$ changes nothing but
  would allow sectoral heterogeneity. **Action:** one parenthetical in the
  surrounding sentence.

# Cut

- **§2 L9 (164 words) --- the OECD skills argument.** The largest single
  footnote in the manuscript, and the only one that carries a full argument
  rather than a qualification: green occupations are mostly existing
  occupations, therefore the program does not demand a new skill frontier. It
  is also the footnote that defends the single-factor design against the
  competence objection, which is why the author may want it in the body
  instead. **Action:** promote its opening sentence into the §2 paragraph and
  cut the rest, or leave it whole. Do not trim it to half an argument.
- **§3 L59 (31 words).** States that under demand-only shocks the elastic model
  collapses to the conventional CGE in its single-wage form. The same fact is
  stated in the §3 L36 footnote and proved in Appendix A. **Action:** delete;
  the two remaining statements already carry it.
- **§6 L20 (28 words).** Repeats four numbers that Appendix E prints: the
  $0.20$ and $0.16$ fall under preference reallocation and tax finance, the
  $+0.12$ rise under external finance, and the standard CGE's $0.16$. This is
  exactly the case rule 5 describes --- estimation output restated in running
  text. **Action:** delete the numbers, keep the pointer to Appendix E.

# Keep

Nine footnotes do work no sentence can do. Each states a definition, a
convention or a distinction the tables depend on, and each is short enough not
to read as displaced prose.

- **§2 L7 (25 words).** Why the inter-industry table is extended with wages and
  consumption shares; the method reference for the whole paper.
- **§3 L36 (65 words).** The single-wage form of the elastic closure and why the
  sectoral form is the one used --- it pre-empts the objection that the
  modelled variety is a choice of exposition.
- **§4 L68 (29 words).** Composition of the Törnqvist price index used for
  deflation; the tables are not checkable without it.
- **§4 L74 (30 words).** When a fixed-weight index coincides with the reported
  measure; this is the property the measure discussion turns on.
- **§4 L78 (79 words).** Public valuation at $\lambda^{G}$ and the stance that
  $\lambda^{G}=1$. A stance stated in the body would read heavier than in a
  footnote, and the welfare table depends on it.
- **§5 L15 (47 words).** Why the elastic column is reported as an interval
  rather than at one elasticity; the bracket convention has no other home.
- **§6 L13 (43 words).** Value-added versus gross output in the
  Kaldor--Verdoorn bracket, with the sector where it bites.
- **§6 L47 (73 words).** The augmented household budget with the retained profit
  share; a definition.
- **(B) L226 (18 words).** $F$ clears the external account while F1--F3 name
  financing closures. Cheap, and it prevents a real reading error.

# Ruled out by the author

Recorded so the next pass does not raise them again.

- **The decimal precision in Appendix C** (three and four places, e.g.
  $0.683$, $0.1199$). The author's verdict: this is a data appendix and the
  precision is appropriate there. Not a shortening candidate.
- **Spelled-out percentages.** "Five percent" and "two percent" in §5 and the
  caption unit labels ("all in percent") stay as written.
- **The number density of Appendix C.** Accepted; the appendix exists to make
  the calibration auditable.

# Revision log

- **v1 (2026-10-24).** First inventory: 19 footnotes, 1,068 words, ten
  candidates at 659 words. Taken from `cf53520` after the mechanical cleaning
  pass; the ruled-out items record the author's decisions of the same date.
