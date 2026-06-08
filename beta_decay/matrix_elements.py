"""
Transition matrix elements for beta decay.

Computes Gamow-Teller (GT) and first-forbidden (FF) matrix elements
from one-body density matrix (OBDM) elements and single-particle matrix elements.

Reference: Behrens & Bühring, "Electron Radial Wave Functions and Nuclear
           Beta-Decay" (1982); Walecka (1975).
"""

import numpy as np
from .wave_function import NuclearState, SingleParticleState


def sp_gt_reduced_me(p: SingleParticleState, n: SingleParticleState) -> float:
    """
    Single-particle reduced Gamow-Teller matrix element:
        <jp || sigma || jn>  (same l required)

    Uses standard shell-model formula:
        <j=l+1/2 || sigma || j=l+1/2> =  sqrt(2(2j+1))
        <j=l+1/2 || sigma || j=l-1/2> = -sqrt(2j+1)  (j = l+1/2)
        <j=l-1/2 || sigma || j=l+1/2> = -sqrt(2j'+1) (j'= l+1/2)
        <j=l-1/2 || sigma || j=l-1/2> = -sqrt(2(2j+1))
    """
    if p.l != n.l:
        return 0.0

    l = p.l
    jp2 = p.j2
    jn2 = n.j2

    if jp2 == 2 * l + 1 and jn2 == 2 * l + 1:
        return np.sqrt(2 * (jp2 + 1))
    elif jp2 == 2 * l + 1 and jn2 == 2 * l - 1:
        return -np.sqrt(jp2 + 1)
    elif jp2 == 2 * l - 1 and jn2 == 2 * l + 1:
        return -np.sqrt(jn2 + 1)
    elif jp2 == 2 * l - 1 and jn2 == 2 * l - 1:
        return -np.sqrt(2 * (jp2 + 1))
    return 0.0


def gamow_teller_strength(state: NuclearState, Ji2: int) -> float:
    """
    Gamow-Teller transition strength:
        B(GT) = |M_GT|^2 / (2*Ji + 1)

    M_GT = sum_{pn} <p||sigma||n> * rho_{pn}

    Parameters
    ----------
    state : NuclearState
        Final nuclear state containing OBDM elements.
    Ji2 : int
        2 * spin of initial state.

    Returns
    -------
    float
        B(GT) value (dimensionless).
    """
    M_GT = 0.0
    for (p_st, n_st), rho in state.items():
        M_GT += sp_gt_reduced_me(p_st, n_st) * rho
    return M_GT ** 2 / (Ji2 + 1)


def gt_matrix_element(state: NuclearState) -> float:
    """Return raw GT matrix element M_GT (before squaring)."""
    M_GT = 0.0
    for (p_st, n_st), rho in state.items():
        M_GT += sp_gt_reduced_me(p_st, n_st) * rho
    return M_GT


# ---------- First-forbidden matrix elements ----------

def sp_ff_r_me(p: SingleParticleState, n: SingleParticleState,
               r_mean: float) -> float:
    """
    Single-particle first-forbidden radial matrix element <r>.
    Approximated as r_mean * delta_{lp, ln±1} * angular_factor.

    r_mean: mean radial matrix element (in fm), typically R_0 * A^{1/3}
    Parity must change: (-1)^(lp+ln) = -1  => lp + ln = odd.
    """
    if abs(p.l - n.l) != 1:
        return 0.0
    # Simple estimate: <r> ≈ r_mean for allowed l-combinations
    return r_mean


def first_forbidden_xi(state: NuclearState, r_mean: float) -> float:
    """
    First-forbidden matrix element xi (rank-0 part):
        xi = sum_{pn} <p||r Y1 sigma||n> * rho_{pn}

    Approximated using <r> ~ r_mean and angular recoupling.
    """
    xi = 0.0
    for (p_st, n_st), rho in state.items():
        xi += sp_ff_r_me(p_st, n_st, r_mean) * rho
    return xi


def first_forbidden_zeta(state: NuclearState, r_mean: float) -> float:
    """
    First-forbidden matrix element zeta (rank-1 part):
        zeta = sum_{pn} <p||r Y1||n> * rho_{pn}
    """
    zeta = 0.0
    for (p_st, n_st), rho in state.items():
        zeta += sp_ff_r_me(p_st, n_st, r_mean) * rho
    return zeta
