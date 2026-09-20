---
title: "Grouping-rule evidence: capacity signals for the sectoral rigidity set"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-20"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [grouping-rule, evidence, vacancies, overtime, backlogs, preregistration, labour]
last-updated: September 2026
---

**Version 1** (September 2026)

Evidence annexe to the grouping-rule decision point of
`docs/WORKPLAN_SENSITIVE_PRICES.md` (v7) and the ADR-0022 open item ("a rule
is preferred to a hand-picked list, to forestall the tuning objection").
The executed `matrix_5x3_v9` two-group cells carry two *rules* (the seven
programme sectors; the largest half by baseline employment) whose choice
flips the sign of the employment effect. This note collects the external,
pre-run evidence that a *principled* third rule can be stated from, and
states that rule. Nothing here has been run; the rule is for ratification
and then preregistration with the design.

# The rule (candidate, for ratification)

> **A sector's labour supply is rigid (`eta_s,i = 0`) iff (i) its programme
> incidence share is at least 1 % of the programme total AND (ii) its
> skilled-worker-shortage signal reaches the goods-producing average or
> above.**

Rationale, in one sentence: the rigidity of a sectoral labour supply for
the programme's purposes is the *elasticity* of its labour market, and the
best available external proxies for short-run inelasticity are persistent
vacancies and reported skills shortages — not the level of capacity
utilisation, which is cyclical and was falling in 2024.

The rule is computable before any run: (i) is read off the 2024 impulse
data, (ii) off the 2024 shortage reports below. Its output — the rigid set
— is a model *input*, not an artefact of the results.

# Programme incidence (i) — from the data

From `cbase2/data_raw/impulses.csv`, 2024 row (71 VGR sectors, G0 =
40,300 EUR m): the programme's incidence is concentrated in seven sectors.

| Sector | 2024 impulse share |
| --- | ---: |
| Specialised construction works | 63.4 % |
| Rubber and plastics products | 9.4 % |
| Ceramic products, processed stone and clay | 8.0 % |
| Chemicals and chemical products | 7.6 % |
| Glass and glassware | 6.1 % |
| Machinery | 4.3 % |
| Electrical equipment | 1.3 % |

Condition (i) selects exactly these seven (the next sector in the ranking
has zero incidence).

# Capacity signals (ii) — 2024 evidence

| Indicator | Value (2024) | Coverage | Source |
| --- | --- | --- | --- |
| Skilled-worker shortage, Bauhauptgewerbe (Q4/2024 survey, October) | 28.9 % of firms report production restricted by skills shortage (27.7 % in July; peak 37.0 % Q4/2022) | WZ broadly building construction + civil engineering + demolition + specialised works (Destatis Bauhauptgewerbe grouping) | KfW-ifo-Fachkräftebarometer December 2024 (KfW Research / ifo Institut), p.2 |
| Skilled-worker shortage, Verarbeitendes Gewerbe aggregate (Q4/2024) | 20.6 % (peak 44.5 % Q3/2022; long-run mean 9.7 %) | manufacturing broad | KfW-ifo-Fachkräftebarometer, December 2024 |
| Skilled-worker shortage, Gummi- und Kunststoffwaren (Q4/2024) | 33.1 % | one of the six material sectors (rubber/plastics) | KfW-ifo-Fachkräftebarometer, December 2024 |
| Skilled-worker shortage, Metallerzeugung/-bearbeitung (Q4/2024) | 10.3 % | comparator: metal processing not programme-relevant | KfW-ifo-Fachkräftebarometer, December 2024 |
| Open positions, Baugewerbe (Q4/2024) | 108,000 (national: 1.40 m; Q4/2022 peak ~2.0 m) | construction sector | IAB-Stellenerhebung, IAB-Monitor Arbeitskräftebedarf 4/2024 (IAB-Forum, 14 Mar 2025) |
| National Vakanzrate (Q4/2024) | 3.2 % (immediately open positions per 100 demanded) | all sectors; 2.9 % in Q3/2024 | IAB-Monitor Arbeitskräftebedarf 4/2024 |
| Order backlog, Bauhauptgewerbe (June 2024) | 3.7 months overall (Hochbau 3.5, Tiefbau 4.1, sonstiger Tiefbau 4.5) | ifo Konjunkturumfragen | ifo / bauletter 1 Jul 2024; declining from 4.8 months (Jan 2022) |
| Working-time accounts, Baugewerbe (2023 -> 2024) | 20 % -> 62 % of firms (Hesse), the largest rise of any sector | IAB-Betriebspanel 2024 | IAB-Betriebspanel Report Hessen 2024 (HNA/iwak) |
| Personnel demand persistence, Baugewerbe (H1 2024) | 54 % of firms with hiring need (3rd highest); of those, 72 % left positions unfilled | IAB-Betriebspanel, Brandenburg 2024 | Betriebspanel Brandenburg 2024 (IAB) |

# What the indicators mean, in plain words

The four sources are three national series and one pair of regional panels.
None of them measures a labour-supply elasticity; together they triangulate
*how hard it is to get workers quickly* in the programme's sectors, which
is the short-run margin the rigidity assumption (`eta_s,i = 0`) is about.

- **IAB-Stellenerhebung** (Institute for Employment Research, quarterly
  establishment survey): firms report all positions they are trying to
  fill, including the ones not registered at the employment agency. Two
  quantities matter: the *stock* of open positions (Baugewerbe 108,000 in
  Q4/2024, of 1.40 million nationally) and the **Vakanzrate** — immediately
  open positions per 100 workers demanded (employed plus open), the
  standard national tightness measure (3.2 % in Q4/2024). Reading: even
  after a year in which vacancies fell 19 % (the 2024 downturn),
  construction still holds a top-five stock of positions that firms keep
  failing to fill — a *supply* symptom, not a demand symptom.
- **KfW-ifo-Fachkräftebarometer** (KfW Research / ifo Institut, quarterly
  survey of about 9,000 firms across manufacturing, the building main
  trades, trade and services): one question, "is your business currently
  restricted by a lack of skilled workers?"; the number is the **share of
  firms** answering yes. Q4/2024: 28.9 % in the Bauhauptgewerbe against
  20.6 % in manufacturing broadly (plastics — rubber/plastics, one of the
  six material sectors — 33.1 %, the most affected manufacturing branch in
  the report). Construction runs at roughly 1.4x the manufacturing rate at
  a moment when manufacturing itself is slack. Note the measure is
  firms' self-reported current restriction and needs a cyclical caveat
  (it fell from a 37.0 % peak in Q4/2022).
- **ifo Konjunkturumfragen — Auftragsbestand** (order backlog in months):
  firms state how many months of work their order books cover at current
  capacity. 3.7 months (June 2024) means: were new orders to stop,
  construction firms still have 3.7 months of scheduled work. This is the
  demand-side/capacity-cushion indicator, and it fell from 4.8 months
  (early 2022). It is in the table deliberately, as the *counter-signal*:
  it documents that the 2024 story is not "order books bursting", so the
  rigidity claim must rest on the labour-elasticity signals (the other
  three), not on utilisation.
- **IAB-Betriebspanel** (annual establishment survey): two regional 2024
  state reports (Hesse, Brandenburg) are used as directional evidence, not
  national estimates. (i) Working-time accounts in the Baugewerbe jumped
  from 20 % to 62 % of firms 2023->2024 in Hesse (37 % -> 58 % West): when
  order books are volatile, firms keep headcount stable and buffer through
  *hours* — hours are the flexible margin, employment the rigid one, which
  is precisely the vertical-supply assumption. (ii) In Brandenburg, 54 % of
  construction firms had a hiring need in H1/2024 (third-highest of all
  sectors) and 72 % of those left positions unfilled — recruitment
  failures persist even in a weak year.

# What the rule does with them

The model's rigidity is a property of the *slope* of a sector's labour
supply (its elasticity), never measured by any of these series. The rule
therefore only uses them to *select the set*: a sector is rigid when it
concentrates the programme (incidence >= 1 %) **and** sits at or above the
goods-producing average on the skills-shortage signal. The value of
`eta_s` itself comes from the supply arm (ADR-0023), not from this table.
What the table establishes, before any model run, is that the seven
programme sectors — not some hand-picked list — are the sectors whose
labour markets show the shortage pattern the rigidity assumption
formalises.

# What the table says, and its limits

- **The capacity-pressure reading survives the 2024 downturn — with a
  caveat.** Order backlogs fell from the 2022 peak (4.8 -> 3.7 months) and
  open positions fell 19 % y/y, i.e. the 2024 demand side was weak. But the
  signals the rule uses are the *elasticity* signals: skills-shortage
  reports and vacancy persistence. Construction still reports 28.9 %
  (8.3 pp above manufacturing), plastics 33.1 %, and 108,000 open
  construction positions; working-time accounts doubled in a year — firms
  buffer volatile order books through hours, which is consistent with
  headcount that cannot expand at will.
- **Provisional output of the rule.** With the collected signals, condition
  (ii) is firmly met by specialised construction and rubber/plastics (both
  above the manufacturing average), plausibly met by ceramics, chemicals
  and glass, and *unverified* for machinery and electrical equipment (2024
  sub-branch shortage shares not retrieved). A defensible pre-run
  specification registers **two** evidence-rule variants: the five-sector
  core {specialised construction, plastics, ceramics, chemicals, glass}
  (94.5 % of programme incidence) and the full seven under a weaker
  threshold (>= 20 %, i.e. the manufacturing average). The five-sector set
  is the default; the seven-sector set tests the threshold's sensitivity.
- **Mapping caveat.** The model's sectors are VGR (Destatis 2019 IO
  classification); the reports above use WZ/NACE groupings. "Specialised
  construction works" ~ WZ 43, and the Destatis Bauhauptgewerbe indicator
  used by ifo/KfW covers WZ 41.2, 42.x, 43.1, 43.9 — a superset. The
  material suppliers map to their WZ classes up to aggregation. This is
  approximation at the 2-digit level, recorded rather than hidden.
- **Vintage.** All signals are 2024, matching the programme vintage
  (impulses.csv 2024 row, G0 = 40,300 EUR m at 2019 prices). The shortage
  shares are firm-level survey shares (KfW-ifo), not labour-market
  elasticities; they are the best public proxies and are used only to
  *select* the rigid set, never as estimates of `eta_s` itself.
- **Not a claim about wage policy.** As in ADR-0022: a welfare difference
  across rules is not a statement about who should earn what; the rule only
  answers which sectors' supply curves are treated as inelastic.

# Next steps

1. Ratification of the rule wording (and the five-vs-seven threshold).
2. The grouping-rule batch design (rules: uniform / programme / evidence-5 /
   evidence-7 / largest-half, x F1/F3 at the settled `eta_s` level), rows
   `planned` in `registry/scenarios.csv`, preregistered (ADR-0006) citing
   this note as the rule's evidence, executed, and the employment sign
   compared across rules — the preregistered signature is that the *rule*,
   not the elasticity level, is the sign-relevant choice.
3. The same evidence ranges become the probability measure of the sectoral
   Sobol (IDEA-0003).

# Revision Log

- **Version 1** (September 2026)