---
title: "IDEA-0008: one index convention for the kernel --- i the subject, u an input, j a user"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-10-23"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, notation, kernel, appendix-b, revision]
last-updated: October 2026
---

**Version 1** (October 2026)

Status: `scoped`.

# The question

Whether the kernel can be written under a single index convention, what such a
convention can and cannot bind, and what migrating to it costs. Recorded because
B.1 states the convention in prose --- "$\Omega_{iu}$ the share of input $u$ in
sector $i$'s bundle" --- and then contradicts it six paragraphs later, and because
the contradiction is invisible to a reader who has learned only the prose.

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

Each letter names a role relative to the equation being written, never a sector's
identity.

| Letter | Role | Positions |
| --- | --- | --- |
| $i$ | the sector the equation is about | fixed index of every block: $P_i$, $c_i$, $p_i$, $y_i$, $L_i$, $\alpha_i$, $A_i$, $c_i^h$, $m_i$, $\mu_i$ |
| $u$ | an input good | summed index inside a bundle; its price $p_u$; the share $\Omega_{iu}$; the supplying slot of $x_{ui}$; final use $z_u$ |
| $j$ | a sector as user, when it is not the subject | summed index of a sum over buyers; the user-side arguments $A_j$, $y_j$, $a_j$ |

What follows mechanically, and therefore need not be restated in the text:

- $\Omega_{iu}$ reads owner first, input second, and its row sums over the second index equal one, $\sum_u \Omega_{iu} = 1$.
- $x_{ui}$ reads supplier first, user second.
- $\Omega_{iu} = p_u x_{ui}/TC_i$ ties the two, so $\Omega$ is the transpose of the column-normalised flow table.
- A user's scale enters as $a_j y_j$: intensity times output.

# Where the rule cannot bind

$i$ cannot be role-pure, and the paper is better off saying so once than
pretending otherwise: in $P_i$ the subject buys, in $\mathit{int}_i$ the subject
is the good demanded, in $y_i$ it sells. In a square model every sector holds
both roles inside the same block. The rule therefore disciplines only the
auxiliary letters $u$ and $j$. That is its honest scope, and it is why the rule is
a *slot* rule rather than a naming scheme: stated as "i is the origin" it would be
false at once.

# The alternative, and its price

The route considered and rejected on 2026-10-23: $u$ = user, $i$ = input, which
would make final use $z_i$ and align the letters with "U = user". It fails in one
place. The price index of the sector whose cost is being computed is $P_i$, and
the same $i$ sits in $p_i = c_i$; if $i$ means input, that sector needs two
letters inside a single equation --- $c_u$ written in terms of $P_u$ with $u = i$.
The duplication then propagates to every $p_i$, $y_i$, $c_i$, $L_i$, $\alpha_i$,
$A_i$ and to the household block, where $\beta_i$ and $c_i^h$ already index goods
by $i$. Cost: all of Appendix B, Appendix H and the Section 4 kernel. Benefit:
none the slot rule does not already deliver.

# The violations that remain

Measured on 2026-10-23 in `revision/sections/(B)model-closures-proofs.tex`. Four
lines hold every case; in each of them $u$ is a buyer, and each is repaired by
relabelling the summed index to $j$ together with its dependent arguments
$p_j$, $A_j$, $P_j$, $a_j$, $y_j$.

| Line | Today | Under the rule |
| --- | --- | --- |
| 315 | $\mathit{int}_i = p_i^{-\theta} \sum_u \Omega_{ui} \dots a_u y_u$ | $\sum_j \Omega_{ji} \dots a_j y_j$ |
| 485 | $\sum_u \Omega_{ui} a_u y_u$ | $\sum_j \Omega_{ji} a_j y_j$ |
| 493 | $\sum_i \Omega_{ui} a_u = a_u$ | $\sum_i \Omega_{ji} a_j = a_j$ |
| 495 | $\alpha_u$ per unit of $y_u$ | $\alpha_j$ per unit of $y_j$ |

# Cost and sequence

Four lines, six token edits, no equation changes, no cross-reference moves, no
numerical consequence. Two orders are defensible: relabel first and rebuild, or
append the mirror sentence to the Notation paragraph and leave the drift, in which
case the reader is warned once and the six sites become tolerable. The first is
cheap enough that the second is only worth choosing if the text is being frozen at
the same time.

# What it buys, and what it does not

It removes the only systematic trap in the kernel: after the migration every
$\Omega$ in the manuscript reads [owner][input] and every $x$ reads
[supplier][user], so a reader who has learned one object has learned all of them.
It does not make the letters carry sector identities, and it does not touch the
ambiguity that no notation can remove --- that the subject $i$ is buyer in one
equation and seller in the next.

# Adjacent findings fixed the same day

- Section 4 defined $\Omega$ transposed, as $\Omega_{ui} = p_u x_{ui}/TC_i$, against nine uses of $\Omega_{iu}$ elsewhere; corrected, with a parenthetical recording that $\Omega$ is the transpose of the column-wise table.
- Section 4's calibration sentence named the measured shares $\Omega_{ui}$; now $\Omega_{iu}$, consistent with its own definition two paragraphs earlier.
- Appendix H wrote $x_{ju}$ and $\Omega_{uj}$ with $u$ as the user and $j$ as the input; aligned to $x_{ui}$ and $\Omega_{iu}$.
- $a_u = A^{\mathrm{bill}}_u/\lambda_u$ carried a symbol, $A^{\mathrm{bill}}$, that occurred nowhere else in the manuscript; replaced by a direct definition, "the value of domestic intermediates used per unit of gross output".
- $g_i = kG_0\psi_i$ carried a factor, $k$, that is one in every matrix cell; now $g_i = G_0\psi_i$, with a $k$-fold program written $kg$ where scaling is invoked.

# Open decisions for the operator

- Apply the six relabellings, or state the mirror once in the Notation paragraph and leave them?
- If applied, does $a$ keep the letter $u$ in its vector definition, where it is a sector-level intensity, while $a_j$ appears inside sums over users? Recommendation: yes, and write the vector bare as $a$ wherever no sum is running.
- Does $z_u$ keep $u$? Recommendation: yes --- it stands beside $y_u$ and $x_{ui}$ in the row identity, where $u$ is the supplying sector.

# Revision Log

- **Version 1** (October 2026) --- recorded the convention, the scope it cannot cover, the rejected alternative, the six remaining sites and the fixes made the same day.
