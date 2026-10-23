---
title: "IDEA-0008: one index convention for the kernel --- i the subject, u an input, j a user"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-10-23"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, notation, kernel, appendix-b, revision]
last-updated: October 2026
---

**Version 2** (October 2026)
**Version 1** (October 2026)

Status: `scoped`.

# The question

Whether the kernel can be written under a single index convention, what such a
convention can and cannot bind, and what migrating to it costs. Recorded because
the closure-systems subsection describes "$\Omega_{iu}$ the share of input $u$
in sector $i$'s bundle", then uses $u$ for a buyer in intermediate demand.
The matrix slots remain consistent, but the letter changes role without warning.

# What drifts today

Two objects disagree about what $u$ means, and each is correct under its own slot
rule:

| Object | Reading | $u$ is |
| --- | --- | --- |
| $P_i(p) = (\sum_u \Omega_{iu} p_u^{1-\theta})^{1/(1-\theta)}$ | the bundle of $i$, summed over its inputs | an input |
| $\mathit{int}_i = p_i^{-\theta} \sum_u \Omega_{ui} p_u^{\epsilon} A_u^{\epsilon-1} P_u^{\theta-\epsilon} a_u y_u$ | demand for good $i$, summed over its buyers | a user |

The second is not an error: $\Omega_{ui}$ carries the owner in first position, so
$u$ is the owner there. What fails is the letter, not the slot. A reader who
carries "u is an input" out of the price index meets $A_u$, $y_u$, $a_u$ in the
demand equation and has to re-derive the rule in the middle of a proof.

# The rule

In production-network expressions, use the following local convention. Letters
name roles within an expression, rather than permanent sector identities.

| Letter | Role | Positions |
| --- | --- | --- |
| $i$ | the sector the equation is about | free index in sector equations such as $P_i$, $c_i$, $p_i$, $y_i$, $L_i$, $\alpha_i$, $A_i$, $c_i^h$, $m_i$, $\mu_i$ |
| $u$ | an input or supplying good | summed index inside a bundle; its price $p_u$; the share $\Omega_{iu}$; the supplying slot of $x_{ui}$; row output $y_u$ and final use $z_u$ |
| $j$ | a sector as user, when it is not the subject | summed index of a sum over buyers; the user-side arguments $A_j$, $y_j$, $a_j$ |

The matrix orientations and calibration normalization must be stated separately;
they do not follow from the choice of letters:

- $\Omega_{iu}$ reads owner first, input second, and its row sums over the second index equal one, $\sum_u \Omega_{iu} = 1$.
- $x_{ui}$ reads supplier first, user second.
- At baseline, $\Omega_{iu} = p_{u,0}x_{ui,0}/TC_{i,0}$, where $TC_{i,0} = \sum_v p_{v,0}x_{vi,0}$ is the domestic intermediate bill, excluding primary inputs, imported intermediates and intermediate-use taxes. Thus $\Omega$ is the transpose of the baseline domestic flow table after normalizing each user column by its domestic bill.
- A user's scale enters as $a_j y_j$: intensity times output.

The weights $\Omega$ stay fixed after calibration. For the CES bundle above,
the current expenditure share is
$s_{iu}(p) = \Omega_{iu}p_u^{1-\theta}/\sum_v\Omega_{iv}p_v^{1-\theta}
= \Omega_{iu}(p_u/P_i)^{1-\theta}$, with the Cobb--Douglas limit at
$\theta=1$. It generally differs from the calibrated weight when relative
prices change. The flow-share identity therefore defines baseline weights,
not an identity that fixes current shares in every equilibrium.

# Where the rule cannot bind

$i$ cannot be role-pure, and the paper is better off saying so once than
pretending otherwise: in $P_i$ the subject buys, in $\mathit{int}_i$ the subject
is the good demanded, in $y_i$ it sells. In a square model every sector holds
both roles inside the same block. The role preference for $u$ and $j$ applies
where a production-network expression distinguishes suppliers from users.
The slot meanings of $\Omega$ and $x$ remain binding regardless of letters.

Generic sector sums and definitions may use ordinary dummy indices. In
particular, the household normalizers $\sum_j\beta_j p_j^{1-\sigma}$ and
$\sum_j\beta_j ds_j p_j^{1-\sigma}$ sum over consumption goods; $j$ does
not denote a buying firm there. Standalone definitions such as $a_u$ may
index sectors generically, provided the local meaning is explicit. Writing
the vector as $a$, or defining its generic component as $a_i$, avoids an
unnecessary switch just before a sum over users.

# The alternative, and its price

An alternative is $u$ = user and $i$ = input, making final use $z_i$ and
aligning the letters with "U = user". This is mathematically valid: one can
write $P_u$ in terms of input prices $p_i$, and $p_u=c_u$, without imposing
$u=i$. Version 1's claim that the alternative requires that equality was
incorrect.

The tradeoff is editorial. Applying that alternative throughout would rename
many existing sector equations and their dependent arguments in Section 4
and Appendices B and H. The proposed local convention preserves the familiar
subject index $i$ and changes fewer expressions. That lower migration cost
is the reason to prefer it; the alternative has no mathematical defect.

# The violations that remain

Inspected on 2026-10-23. These are the production-network expressions to
align; line numbers refer to the source at that date. Rename each index
together with all dependent arguments. This is a scoped inventory, not a
claim that every use of $u$ or $j$ in the manuscript has been audited.

| Source and line | Today | Under the rule |
| --- | --- | --- |
| Appendix B, 291–292 | sector coefficients $a_u$ | define $a_i$ or the vector $a$ to avoid a local role switch |
| Appendix B, 315–316 | $\mathit{int}_i = p_i^{-\theta} \sum_u \Omega_{ui} \dots a_u y_u$ | $\mathit{int}_i = p_i^{-\theta}\sum_j\Omega_{ji}p_j^{\epsilon}A_j^{\epsilon-1}P_j^{\theta-\epsilon}a_jy_j$ |
| Appendix B, 482 | $P_u=1$ in the intermediate-demand proof | $P_j=1$ |
| Appendix B, 485 | $\sum_u \Omega_{ui} a_u y_u$ | $\sum_j \Omega_{ji} a_j y_j$ |
| Appendix B, 493 | $\sum_i \Omega_{ui} a_u = a_u$ | $\sum_i \Omega_{ji} a_j = a_j$ |
| Appendix B, 495 | $\alpha_u$ per unit of $y_u$ | $\alpha_j$ per unit of $y_j$ |
| Appendix H, 23 | purchaser prices $p_j$ in the definition of $P_i$ | input prices $p_u$ |

The source paths are `revision/sections/(B)model-closures-proofs.tex` and
`revision/sections/(H)pricing-door.tex`. The household sums using $j$ are
explicit exceptions, so they require no relabelling.

# Cost and sequence

This is a small manuscript edit spanning several expressions and their
surrounding definitions. A consistent renaming has no mathematical or
numerical consequence and requires no kernel change or new run generation.

If adopted, state the local convention and its generic-index exceptions,
qualify the flow-share definition as baseline calibration, and apply the
inventory above. Then inspect all active manuscript uses of $\Omega$, $x$,
$u$ and $j$ for orientation and locally consistent dependent arguments.
Rebuild the revision with `latexmk -pdf` and inspect the affected equations.
This note specifies that future edit; it does not implement it.

# What it buys, and what it does not

It reduces role switching between bundle-price equations and intermediate
demand. The slot rule remains $\Omega$ = [owner][input] and $x$ =
[supplier][user]; the letter preference helps readers track those roles.
It does not assign permanent sector identities to letters or establish that
all other notation issues have been resolved. The subject $i$ can still buy
in one equation and sell in the next.

# Adjacent findings fixed the same day

Version 1 recorded the following adjacent edits. The baseline-share
qualification proposed above remains to be applied to the manuscript.

- Section 4 defined $\Omega$ transposed, as $\Omega_{ui} = p_u x_{ui}/TC_i$, against nine uses of $\Omega_{iu}$ elsewhere; corrected, with a parenthetical recording that $\Omega$ is the transpose of the column-wise table.
- Section 4's calibration sentence named the measured shares $\Omega_{ui}$; now $\Omega_{iu}$, consistent with its own definition two paragraphs earlier.
- Appendix H wrote $x_{ju}$ and $\Omega_{uj}$ with $u$ as the user and $j$ as the input; aligned to $x_{ui}$ and $\Omega_{iu}$.
- $a_u = A^{\mathrm{bill}}_u/\lambda_u$ carried a symbol, $A^{\mathrm{bill}}$, that occurred nowhere else in the manuscript; replaced by a direct definition, "the value of domestic intermediates used per unit of gross output".
- $g_i = kG_0\psi_i$ carried a factor, $k$, that is one in every matrix cell; now $g_i = G_0\psi_i$, with a $k$-fold program written $kg$ where scaling is invoked.

# Open decisions for the operator

- Adopt the local convention and apply the inventory above? Recommendation: yes, together with the baseline-share qualification and a manuscript-wide orientation check.
- How should the standalone domestic-bill coefficients be introduced? Recommendation: write the vector as $a$ and define $a_i$ for a generic sector; use $a_j$ inside sums over users. Generic $a_u$ is permitted but less clear in the adjacent definition.
- Does $z_u$ keep $u$? Recommendation: yes --- it stands beside $y_u$ and $x_{ui}$ in the row identity, where $u$ is the supplying sector.

# Revision Log

- **Version 2** (October 2026) --- narrowed the letter convention to local production-network roles, stated generic-index exceptions, qualified the flow-share identity at baseline, corrected the argument against the alternative, expanded the migration inventory, and specified manuscript verification without claiming implementation.
- **Version 1** (October 2026) --- recorded the convention, the scope it cannot cover, the rejected alternative, the six remaining sites and the fixes made the same day.
