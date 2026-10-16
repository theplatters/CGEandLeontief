#!/usr/bin/env python3
"""Does a preference-consistent re-weighting rescue W_b under F1?

The proposal: if F1 is a preference shift, re-base the bundle aggregator on the
shifted preferences -- then it should reproduce the exact money-metric W_a'.

It does not. A preference change destroys the cardinal scale, so there is no
canonical ratio across the two regimes. This fixture quantifies the spread.

    python3 welfare-preference-shift.py
"""

import numpy as np

SIGMA = 0.9
SAVE = 0.0
OMEGA = np.array([0.6, 0.4])          # baseline consumption shares
TILT = np.array([2.0, 0.5])           # the F1 preference shifter d_i
E0, E1 = 1.0, 1.05                    # income effect is +5.0 %


def ces(q, q0, w):
    """CES quantity index with weights w (q0-relative), kernel convention."""
    r = q / q0
    return (w * r ** ((SIGMA - 1) / SIGMA)).sum() ** (SIGMA / (SIGMA - 1))


def utility(c, wts):
    """Direct CES utility, weights normalised to sum to one.

    With unit prices this is exactly real income E -- the money-metric
    normalisation -- which is why this reading collapses onto E1/E0.
    """
    return (wts ** (1 / SIGMA) * c ** ((SIGMA - 1) / SIGMA)).sum() ** (SIGMA / (SIGMA - 1))


def main():
    c0 = (1 - SAVE) * OMEGA * E0              # baseline demands, p = 1
    alpha = OMEGA * TILT
    alpha = alpha / alpha.sum()               # shifted preference weights
    c1 = (1 - SAVE) * alpha * E1              # F1 cell demands, p = 1

    exact = E1 / E0
    w_quant = c0 / c0.sum()                   # what ces_quantity_index uses
    w_pref = OMEGA * TILT ** (1 / SIGMA)      # preference-consistent weights
    w_pref = w_pref / w_pref.sum()

    print(f"c0 = {c0},  c1 = {c1},  alpha = {alpha}")
    print(f"real-income effect  E1/E0             = {exact:.6f}   ({100*(exact-1):+.2f} %)")
    print(f"CES, baseline-quantity weights (W_b)  = {ces(c1,c0,w_quant):.6f}   "
          f"({100*(ces(c1,c0,w_quant)-1):+.2f} %)")
    print(f"CES, preference-consistent weights    = {ces(c1,c0,w_pref):.6f}   "
          f"({100*(ces(c1,c0,w_pref)-1):+.2f} %)")

    # third reading: each state under its OWN preferences, in money-metric
    # units (weights normalised to sum to one, so U = E at unit prices).
    # This collapses identically onto the real-income effect -- a CES
    # household's utility IS its real income.
    u0 = utility(c0, OMEGA)
    u1 = utility(c1, alpha)
    print(f"\nU_omega(c0) = {u0:.12f}   U_alpha(c1) = {u1:.12f}")
    print(f"each state by its OWN preferences      = {u1/u0:.12f}   "
          f"({100*(u1/u0-1):+.2f} %)  <- IS the real-income effect, identically")

    # control: with d = 1 the basket grows proportionally and every
    # homogeneous index agrees with the real-income effect (cf. the F2/F3 rows)
    a = OMEGA / OMEGA.sum()
    c1p = (1 - SAVE) * a * E1
    print(f"control d=1: baseline-weighted        = {ces(c1p,c0,w_quant):.6f}   "
          f"({100*(ces(c1p,c0,w_quant)-1):+.2f} %)  <- matches, as it must")


if __name__ == "__main__":
    main()
