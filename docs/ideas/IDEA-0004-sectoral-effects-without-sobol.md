---
title: "IDEA-0004: sectoral effects without a Sobol, and reviving the paper's sectoral figures"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-20"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, sectoral, incidence, divisia, affected-vs-unaffected, visualization, figures, follow-up]
last-updated: September 2026
---

**Version 1** (September 2026)

Status: `scoped` (the measures, the data basis and the cost are specified; nothing measured).

# The idea

A sectoral analysis of the executed matrix that needs **no** Sobol. The sectoral decomposition of `IDEA-0003` is the expensive route; most of the question it addresses --- *what happens to the 71 sectors across the closures, and which sectors carry the headline* --- can be answered descriptively from artifacts that already exist, at a fraction of the cost, and it should be executed **before** any Sobol design is preregistered, because its output tells the operator whether the expensive block is worth its solves.

The companion point is a figure question. The original submission already had a sectoral visual grammar; the review below asks whether reviving it is worth the re-generation. The answer is yes for two of its exhibits and no for a third.

# What is available now (the data basis)

| Item | State |
| --- | --- |
| Sectoral price and quantity vectors | in `runs/<run_id>/solution.csv`, columns `sector,price,quantity`, 71 rows + header |
| Coverage | all 33 `matrix_5x3_v9` cells have a `solution.csv` on this working copy; **none** of the 33 `matrix_5x3_v10` dirs does (manifests only) |
| Which generation to read | the C1 defect was a *measurement* defect: `gdp_components` collapsed the sectoral wage vector, so the solve is untouched (prices, quantities, employment, consumption identical, ADR-0022 amendment). Sectoral `p` and `y` from the v9 solutions are therefore usable; any **gdp-family** number must be cited to v10 |
| Sectoral wage and employment vectors | **not stored** --- `solution.csv` has no wage or employment column |
| Household consumption by sector | **not stored** --- recoverable only by re-evaluating the demand hook |
| `src/` change | none required for any measure below |

The wage gap matters for exactly one of the five measures and is a small harness change, noted there.

# The measures

## Contribution decomposition (Divisia/Hulten weights)

The first thing to do, and the one that connects the sectoral detail to the published headline exactly rather than approximately.

- **What it answers:** which sectors *contribute* to the aggregate change, as distinct from which sectors *move*. Weight each sector's log-change by its baseline share so that the aggregate is reproduced: `d log X = sum_i omega_i d log x_i`. The weights are the model's own baseline shares, not an assumption.
- **Why it is the right first exhibit:** it is the paper's own subject matter. The whole BeyondHulten argument is that aggregation is non-linear, so the weights --- and hence the sectoral composition --- determine the aggregate. A contribution table makes that visible per cell instead of asserting it.
- **Cost:** reading existing artifacts; no solves.
- **Standard:** this is the standard sectoral reporting convention in CGE/IO work; nothing to defend.

## Affected versus unaffected sectors, graded

- **What it answers:** whether the programme's effect stays inside its seven incidence sectors or leaks into the rest of the economy.
- **The honest form is graded, not binary.** The binary split is informative but its boundary is arbitrary --- it is the grouping-rule problem ADR-0022 leaves open --- and it hides roundabout spillovers by construction. Report the binary headline, then the graded version as the robustness check: `d log x_i` against baseline programme-incidence share, and against import margin `m_i` and labour share.
- **What to expect:** the submission found the opposite signs in the two groups --- the Leontief model expands almost all sectors, the CGE reduces gross output in most *unaffected* sectors. That asymmetry, not the aggregate, is the sectoral finding.
- **Cost:** reading existing artifacts.

## Direct versus indirect (network) split

- **What it answers:** how much of each sector's change is its own programme bundle and how much is input-output roundabout.
- **Why it is informative:** it is the theoretically grounded version of "affected versus unaffected" and it isolates whether the action is the direct bundle or the propagation --- which is the Baqaee-Farhi channel the model is built on.
- **Cost:** needs the A-bill matrix alongside `y`; cheap, no solves.

## Dispersion and cross-cell similarity

- **What it answers:** how uneven the reallocation is, and whether two closures are telling the same sectoral story.
- **Measures:** CV, Gini or max-over-median of `d log x_i` for the dispersion; Spearman rank correlation or cosine similarity of the `d log p` and `d log y` **vectors** between two cells for the similarity. With only five labour rows and three financing columns, a 15-cell similarity matrix is a legible object.
- **Cost:** reading existing artifacts.

## The wage-to-price map

- **What it answers:** how the sectoral wage vector becomes the sectoral price vector --- the pass-through that no prior closure in this model could produce.
- **Why it is the most publishable sectoral exhibit available:** it is the visual form of the labour door ADR-0022 opens. A scatter of `d log p_i` against `d log w_i` across the `eta_s` ladder and the three financing columns shows the mechanism directly, and F1 being distinct from F2 = F3 at every rung is the demand-sensitivity signature.
- **The gap:** the wage vector is not in `solution.csv`. Fixing it is a small `experiments/run.jl` change (write a `wage` column, and an `employment` column while at it) plus a re-run; no `src/` change, so no kernel provenance invalidation. Do this **once**, together with the figure re-generation, rather than as two batches.

# The submission's sectoral analysis, and whether to return to it

## What the paper actually did

Five sectoral exhibits, with the readings it drew from them.

| Figure label | File | What it shows |
| --- | --- | --- |
| `fig:elasticity-gradient` | `pictures/impulse.png` | four-panel one-at-a-time elasticity sweep of real GDP, with Leontief and Cobb-Douglas benchmarks |
| `fig:comparison-sec4` | `pictures/panel.png` | two panels: prices versus consumption for the seven shocked sectors, and quantities versus prices for all sectors --- shocked in red, the remaining in pink, with error bars marking the range over two extreme elasticity vectors |
| `fig:sectoral-changes` | `pictures/diff_lambda_imp.png` | sectoral gross-output changes, shocked sectors unshaded, Leontief against CGE |
| `fig:elasticity-gradient2`, `fig:panel_ls`, `fig:sectoral-changes_ls` | `impulse_ls.png`, `panel_ls.png`, `diff_lambda_imp_ls.png` | the same three under the legacy labour-slack closure |
| `comparison_between_labor_slacks` | `comparison_between_labor_slacks.png` | the slack variants against each other |

The readings the paper drew, which are the things a revival would have to re-establish or overturn:

- **A positive price-quantity correlation** in the CGE response (quadrant I dominant) --- read as the demand-shock signature, prices of high-demand goods rising.
- **Five named second-quadrant outliers** (rising price, falling quantity): *Lumber and Wood*, *Metal ores and Mining*, *Coke and refined petroleum products*, and two shocked sectors, *Chemicals and chemical products* and *Electrical equipment*. The paper explains these as a **cost shock** transmitted through intermediate prices.
- **The group asymmetry:** the Leontief model expands almost all sectors; the CGE reduces gross output in most unaffected sectors. These offsetting movements, not the aggregate, produce the near-zero real-GDP response.
- **Two idiosyncratic sectors:** 33 *building construction works* and 34 *civil engineering works* expand intermediate production at the expense of final output.

## Merit in returning to the visualization

**Verdict: high merit for the price-quantity scatter and the sectoral output chart; no merit for the elasticity-gradient panel; and the return is a re-generation, not a reprint.** Confidence: high on the first two, moderate on the third (it depends on whether the sectoral story becomes a headline dimension of the revision --- ADR-0022 leaves that open).

| Exhibit | Verdict | Reason |
| --- | --- | --- |
| Price-quantity scatter | **revive, upgrade** | the mechanism reversal is now measurable and it is the paper's thesis |
| Sectoral output chart | **revive** | the group asymmetry is the sectoral finding; the colour grammar already exists |
| Elasticity-gradient panel | **retire** | superseded by the aggregate Sobol, which answers the same question with a variance share |
| The three labour-slack exhibits | **retire with the closure** | the legacy slack callback is dropped; the figures belong to it |

The upgrades, in order of value:

- **The mechanism reversal is the exhibit.** The published scatter established a demand-driven signature under a closure with no wage-vector channel and no financing. Under the sectoral-labour closure prices are cost-pushed by the wage vector: a sector with a rising wage should show a rising price and a falling quantity --- the paper's *outlier* quadrant becomes the *expected* location. Running the same scatter across the five closures therefore separates the demand-driven from the cost-driven regimes on one page, which is precisely the measured map the assessment claims.
- **The financing columns are new information.** The published figures are unfinanced --- the shock is pure manna, and the paper says so. F1/F2/F3 as a third axis is something the submission could not show.
- **The error bars can be upgraded, not merely repeated.** They are a two-point min-max range over extreme elasticity vectors --- the crude ancestor of the Sobol. Keep them as a range band, and where the Sobol exists (block B of `IDEA-0003`) replace the band with a variance share. That converts a 2-point range into a decomposition, which is the honest way to present both.
- **The shock magnitude is comparable.** The paper's impulse is 40.337 bn EUR, about 1.3% of GDP; the current programme is `G0 = 40300` in 2019 prices, which ADR-0024 identifies as the same number. A regenerated figure is directly comparable to the published one, which makes the re-run a replacement rather than a new experiment.
- **The named sectors give the figure a text.** The cost-shock sectors are nameable and testable ex ante: under the current model they should be the sectors with the highest intermediate exposure to programme sectors. A table of predicted against measured location would make the scatter self-diagnosing.

The caveats, which are the reason this is a re-generation:

- **Provenance.** The PNGs in `paper/pictures/` come from the retired pipeline --- legacy labour-slack closure, unfinanced shock, pre-ADR-0010/0012/0019 kernel. They cannot be reprinted or lightly re-plotted; the figures must be regenerated from a current generation.
- **The elasticity-gradient panel should not be revived.** It is one-at-a-time, so it cannot show interactions, and the paper's own reading ("no common monotonic relationship ... sensitivity becomes visible mainly when several elasticities are simultaneously very low") is partly an artefact of that design. The Sobol answers the same question with a variance share. Keep the panel only as the historical motivation for the Sobol.
- **The wage-coloured version needs the wage vector**, which `solution.csv` does not carry. One harness change plus one re-run, as noted above.

## The revived figures: specification

What step 5 of the sequence below actually produces. Two figures, one writer extension.

### Figure R1 --- sectoral price-quantity scatter (revival of `fig:comparison-sec4`)

| Element | Specification |
| --- | --- |
| Panels | small multiples: rows = the price-moving regimes (the `BF` wage-vector row, the `BETA` sectoral `eta_s` ladder, a `GAMMA` tilt row), columns = F1 / F2 / F3 |
| Single-wage rows | one shared control panel: ALPHA, scalar BETA, DELTA and the untilted GAMMA pin have `max abs(p-1)` of order 1e-14, so they are one point cloud, not five figures |
| Axes | `d log p_i` against `d log y_i`, one point per sector |
| Colour | the submission's grammar, reused: the seven incidence sectors in one colour, the remaining 64 in another, and the second-quadrant membership in a third so the cost-shock prediction is visible |
| Range band | replace the two-point min-max error bars with the preregistered `eta_s` ladder range (cheap) or the block-B variance share (expensive); the choice must be stated, not left to the plotting code |
| Labels | the named cost-shock sectors labelled; the remaining sectors as unlabelled points |
| Generation | `matrix_5x3_v10` for any gdp-family number; the `p` and `y` vectors are v9-identical |

### Figure R2 --- sectoral output chart (revival of `fig:sectoral-changes`)

| Element | Specification |
| --- | --- |
| Form | sorted bar chart of `d log y_i` by sector, one panel per closure row, F2 as the reference financing column |
| Colour | incidence sectors unshaded --- the submission's convention --- and the remainder shaded |
| Reference line | the aggregate `d log` at the same cell, so the offsetting movements are visible rather than asserted |
| Comparator | the `DELTA` endpoint bar replaces the old Leontief-versus-CGE contrast: it is the model's own IO corner |

### What each figure must establish, and what would falsify it

- R1 must reproduce or overturn the **sign** of the price-quantity correlation, per regime. If the sectoral rows still show a dominant quadrant I, the cost-push story is wrong and the mechanism-reversal claim must be dropped, not softened.
- R2 must reproduce the **group asymmetry** --- contraction in unaffected sectors --- or show that financing removed it.
- Both must place the submission's **named sectors** where an exposure-based prediction puts them; if they do not, the self-diagnosing table is dropped rather than explained away.

### The writer extension and the batch

- Extend the solution writer **once** with `price, quantity, wage, employment, consumption`. The wage column serves the wage-to-price map; the employment column serves the allocation figures; the consumption column is what the submission's **left panel** (prices versus consumption for the incidence sectors) needs, so reviving it costs nothing extra.
- Re-mint the batch as a new generation (`matrix_5x3_v11`). The writer is harness code, not a kernel file, so no provenance invalidation under ADR-0001/0004 --- but the artifacts change, so the batch is re-minted rather than amended. The metric block must reproduce v10; the provenance block necessarily differs (a new commit).
- Prototype R1 on a single cell (`BETA-F2-etas05`) before the full grid. The critical path first: one cell proves the columns, the axes and the mechanism story, and it costs one solve.

# Recommended sequence

1. Contribution decomposition across the executed cells (existing artifacts, no solves).
2. Affected-versus-unaffected, binary headline plus graded robustness.
3. Direct-versus-indirect split.
4. Cross-cell similarity matrix.
5. Extend the solution writer with `price, quantity, wage, employment, consumption`, re-mint one generation, and regenerate the two figures to the specification above --- at which point the wage-to-price map is available from the same columns.
6. Only then decide whether `IDEA-0003` block B is worth its solves.

Steps 1 to 4 are reading exercises and can be done in one session. Step 5 is a batch. Step 6 is a decision, not a task.

# Anchored in

- `paper/submission-metro/submission.tex` §5 (`subsec:CGE_elasticity`, the inter-paradigmatic comparison) and §6 (`sec:labor`) --- the five sectoral exhibits and their readings.
- `paper/equivalence.tex` Lemma 1 and the sectoral section (why single-wage closures are demand-free in prices, hence why the scatter separates the regimes).
- `docs/decisions/ADR-0022-sectoral-labour-markets.md` (the sectoral family, the C1 amendment, the citation rule for v9 versus v10); `docs/decisions/ADR-0004-runs-are-immutable-manifests.md` (artifacts and citation).
- `docs/grouping_rule_evidence.md` (the seven programme sectors and the incidence shares); `docs/decisions/ADR-0024-programme-total-settled.md` (the shock magnitude).
- `docs/ideas/IDEA-0003-sectoral-sobol.md` (the expensive companion route).

# Revision Log

- **Version 1** (September 2026) --- the measures that need no Sobol, the data basis and its gaps, the submission's sectoral exhibits with a revival verdict, the specification of the two revived figures with their falsifiers and the writer extension, and the recommended sequence.
