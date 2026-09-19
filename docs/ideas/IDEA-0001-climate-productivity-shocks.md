---
title: "IDEA-0001: productivity shocks framed as climate change"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, supply-shock, climate, productivity, identifiability]
last-updated: September 2026
---

**Version 1** (September 2026)

Status: `unexplored`. No specification has been written and nothing has been run.

# The problem it addresses

Every cell of the evaluation matrix is a demand-side shock: the programme
impulse, the F1 preference tilt, the F2 tax and the F3 external financing. There
is no supply-side variation at all, and three coincidences in the current
generation follow from that single fact. BETA reproduces ALPHA to six digits
because with one factor, constant returns and a demand-free price block the
labour-supply elasticity has nothing to bite on; DELTA reproduces GAMMA because
the CES and Leontief technologies differ only once prices move; and no cell
carries a price response at all outside the `eta = 0` row, where the sectoral
wages are solved (measured `max abs(p - 1)` of 0.221 to 0.279 against 2.7e-15 in
the ALPHA and BETA rows).

A supply-side shock is the obvious remedy, and it is exactly what makes the
paper's demand-shock framing awkward: an ad hoc productivity shock bolted onto a
programme-financing exercise reads as a robustness curiosity, and it invites the
question of why it is there.

# The idea

Give the supply side a story the paper already owns. Instead of an unexplained
productivity shock, let the shock be a **climate-change scenario** with sectorally
differentiated incidence, so its incidence vector is argued rather than assumed:

- A **damage** variant: sectoral productivity losses from heat, drought or
  flooding, concentrated in the sectors a physical-impact literature identifies.
- A **transition** variant: sectoral productivity gains in the sectors the green
  programme expands, which is the same subject matter as the programme itself and
  turns the supply side into the second half of the paper's own story rather than
  a foreign perturbation.

The two variants also have opposite signs, which is what makes the exercise
informative: a damage scenario and a transition scenario can bracket the
identifiability question rather than answer it on one side only.

# What it would buy

- Identifiability of `eta_s`: BETA's labour-supply elasticity is currently
  unidentifiable under demand-only shocks, which the assessment records as a
  model-class property rather than a mis-specification.
- A genuine CES-versus-Leontief contrast: DELTA stops being a six-digit copy of
  GAMMA once prices move.
- A price dimension in every labour row, which the matrix currently has only in
  the `eta = 0` row.
- A climate narrative attached to the supply side, consistent with the paper's
  green-investment subject.

# What it would need

- A supply-shock specification in the design schema and in the harness. The
  kernel already accepts one (the first argument of `Shocks`); the harness
  hard-codes a null shock in `build_cell_model`, so no cell has ever carried one.
- An **incidence vector** for the climate shock: which sectors, what magnitude,
  and from where. This is exogenous, with the same epistemic status as the
  programme incidence `psi`, and it must be preregistered.
- Reference data to be collected rather than assumed. Candidate sources for
  sectoral climate damages and for green-transition productivity effects would
  have to be identified with the citation workflow before any calibration; this
  note deliberately names none, because none has been checked.
- A decision on whether the shock enters as a supply-shock multiplier only, or
  whether a second factor is needed for the aggregate effects to be meaningful.
  The larger the intended story, the more the exercise stops being an addition to
  the current model and starts being a different model.

# Risks and objections

- **Framing risk.** A climate shock in a paper about a demand-side programme can
  read as opportunistic. It is defensible only if the incidence vector is argued
  and the scenario is presented as a distinct exercise, not folded into the
  programme's headline numbers.
- **Calibration risk.** Sectoral climate damages are contested and
  heterogeneous; a single incidence vector will not settle that, and the note
  should not be read as claiming that it would.
- **Scope risk.** The cheapest version (a multiplier on `supply_shock`) is a
  robustness exercise. The interesting version (transition-driven productivity
  gains interacting with the programme) is a paper.
- **Timing risk.** This is a revision. Unless a referee asks for a supply side,
  this belongs to the next paper, and the register exists so that it survives the
  revision rather than being rediscovered.

# Where it sits relative to the current work

Door 2 of `docs/VariationinGamma.md` records the mechanical route (a
supply-shock spec in the design and harness) as the substitute for the wage-side
work the paper has decided against. This note is the thematic version of that
route: it does not change the mechanism, it supplies a reason for the shock and a
calibration discipline for its incidence.

ADR-0021's option C is the wage-side complement, kept for the case where the
price response has to be generated from demand-side variation alone; the
workplan in `docs/WORKPLAN_SENSITIVE_PRICES.md` pursues that route instead.

# Revision Log

- **Version 1** (September 2026)