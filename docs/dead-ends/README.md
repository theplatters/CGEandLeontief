# Dead ends

Append-only register of approaches that were tried and abandoned. Its purpose
is to stop agents and humans from re-trying the same failures: every record
states what was tried, how it failed (with evidence), what replaced it, and
the conditions under which it could be revived.

## Rules

- One file per dead end: `DE-####-short-slug.md`, numbered consecutively,
  never reused. The id is referenced from `registry/closures.toml`
  (`dead_ends = [...]`) and from ADRs.
- A record is required when an approach that was taken seriously is
  abandoned — not for typos or ordinary bug fixes. When in doubt, record it.
- Records are append-only: do not rewrite an old record. If a dead end is
  revived, add a new ADR and note in the registry; the DE stays as history.
- A commit that abandons an approach must add the DE record (see `AGENTS.md`).

## Template

```markdown
# DE-#### — Title

- **Status:** dead end — <one-line reason>
- **Recorded:** YYYY-MM-DD
- **Origin:** where/when it was used
- **Related:** ADRs, closure ids, other DEs
- **Evidence:** file:line, diagnostics, run ids, doc sections

## What was tried
## Why it seemed promising
## How it failed
## What replaced it
## Revival conditions
```

## Index

| ID | Title |
| --- | --- |
| DE-0001 | Unnormalized household demand shifter |
| DE-0002 | Dropping a zero-profit equation to close the system |
| DE-0003 | `η = 10⁶` as a fixed-real-wage approximation |
| DE-0004 | Nominal wage in the labour-supply function |
| DE-0005 | Fixed-base quantity sum called real GDP |
| DE-0006 | Variance shares renormalized over main effects / silent incomplete designs |
| DE-0007 | Adopting Baqaee–Farhi (2022) wholesale as the base model |
| DE-0008 | Open-economy calibration without a saving rate (cbase2 v2) |
| DE-0009 | Orphan src files archived |
