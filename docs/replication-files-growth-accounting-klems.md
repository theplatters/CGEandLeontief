---
title: "Replication Files/Growth Accounting_Klems"
covers:
  - "Replication Files/Growth Accounting_Klems/"
reviewed: "cf414f19419f6fa40bc04d0352326807e7a226cc"
---

# Replication Files/Growth Accounting_Klems

This directory is one of the self-contained MATLAB replication arms under
`Replication Files/`. The parent `Read_Me.txt` gives it as the place where the
growth-accounting results of "Section 6.3" are reproduced, alongside the
88-sector GDP-simulation and adjustment-cost arms in sibling folders.[^1] The
whole `Replication Files/` tree, this folder included, entered the repository in
a single `init` commit and has not been edited since; the status board freezes
it as `source-archive` (read-only, 106 files, tree hash `d4b7ed1`), described as
"MATLAB replication code and paper source material for the Baqaee-Farhi (2019)
replication", with `Growth Accounting_Klems/` named explicitly.[^2]

## What the code does

`Growth_Accounting.m` is the directory's only script. It reads the KLEMS
workbook `usa_wk_mar_2017.xlsx` (sheet `KLEMdata`, range `A3:J4422`) and
reshapes the ten columns into 68 sectors by 65 periods, giving gross output,
capital, labour and intermediate use twice: a nominal series used for the
shares and a constant-price series used for the growth rates. From the nominal
levels it forms labour, capital and intermediate shares and a nominal GDP as the
sum of capital and labour income, then computes per-sector Solow residuals two
ways — a full `sTFP` that subtracts all three factor contributions and a
value-added `sVATFP` that subtracts only capital and labour (the intermediate
term is commented out).[^3]

The script then aggregates. Domar weights are gross output over GDP and the
nonlinear aggregate TFP is the Domar-weighted sum of the sectoral residuals;
two approximations serve as checks — a first-order one holding Domar weights at
their base-year values and a second-order one using the Törnqvist average of the
first and last period. The resulting figure plots the cumulated nonlinear,
first- and second-order TFP on a log scale from 1947 to 2013 and is exported to
`Cost_Disease.tex` through `matlab2tikz`; two further in-script figures compare
average growth rates and the per-step Hulten-approximation error.[^4]

## How it is organised and why

The folder is deliberately minimal: one script, one data workbook, the generated
TikZ figure `Cost_Disease.tex`, and three bitmap exports (`Base_year_1947.png`,
`Growth_Rate.png`, `Error.jpg`). It carries no helper functions of its own — the
plotting relies on `matlab2tikz.m`, which the parent `Read_Me.txt` lists among
the utilities that have to be placed on the MATLAB root path, and
`Growth_Accounting.m` is the single entry point that `Read_Me.txt` gives for the
Section 6.3 results.[^5] It is written as a one-shot script rather than a
library because it recreates one paper figure and one set of numbers; the freeze
board keeps the zone read-only "for reference only".[^6]

The script carries the marks of a research working copy: `cd` and `xlsread` are
hard-coded to the original author's Dropbox directory, so it cannot run
unchanged outside that machine.[^7]

## Relationship to the current codebase

The live model is no longer MATLAB. The repository was split: the kernel moved
to the separate OpenCGE Julia package, and this repository keeps the paper, its
run manifests, and the frozen legacy material.[^8] This directory contributes
nothing to `registry/`, no run manifest and no test. The status board records
that it is cited by `bf_replication/REPLICATION_WORKPLAN.md` as the location of
the original MATLAB code, and the workplan maps `Growth Accounting_Klems/` to the
Section 6.3 growth accounting while stating that the Baqaee-Farhi (2020) COVID
replication source is not present in this repository — so the folder is kept as
the source record of that exercise, not as runnable pipeline code.[^9]

[^1]: git:c544907:Replication Files/Read_Me.txt — lists `Growth Accounting_Klems\Growth_Accounting.m` for the Section 6.3 results next to the 88-sector and adjustment-cost runs.
[^2]: status.md — freeze board row `Replication Files` (source-archive, read-only, 106 files, tree hash `d4b7ed1`).
[^3]: git:c544907:Replication Files/Growth Accounting_Klems/Growth_Accounting.m — the `xlsread` range and reshape, the nominal/constant-price columns, and the `sTFP`/`sVATFP` definitions.
[^4]: git:c544907:Replication Files/Growth Accounting_Klems/Cost_Disease.tex — the `matlab2tikz` output with the Nonlinear / First Order / Second Order curves over 1947–2013; the in-script `plot` calls define the average-growth-rate and error figures.
[^5]: git:c544907:Replication Files/Read_Me.txt — names `matlab2tikz.m` among the root-path utilities and `Growth_Accounting.m` as the Section 6.3 entry point.
[^6]: status.md — the freeze board keeps the zone read-only "for reference only".
[^7]: git:c544907:Replication Files/Growth Accounting_Klems/Growth_Accounting.m — `cd` and `xlsread` are hard-coded to `C:\Users\baqaee\Dropbox\Work\Hultens Theorem\Calibration\Growth Accounting_Klems`.
[^8]: decisions/ADR-0028-split-paper-and-opencge-kernel.md — the kernel/model split into the OpenCGE package and this paper repository.
[^9]: status.md — freeze board notes: "Referenced by `bf_replication/REPLICATION_WORKPLAN.md` as the location of the original MATLAB code", and the Baqaee-Farhi (2020) replication source is absent.
