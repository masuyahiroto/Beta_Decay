"""
Quenching factor for Gamow-Teller and first-forbidden transitions.

The free-nucleon GT operator is quenched in nuclei due to short-range
correlations and delta-isobar effects. A standard effective value is:
    g_A^eff = q * g_A,  q ≈ 0.74  (mass-number dependent)

Reference: Wildenthal et al. (1983); Martinez-Pinedo et al. (1996).
"""

from .constants import G_A

# Standard quenching factor for sd/pf shell nuclei
Q_STANDARD = 0.74

# For heavy nuclei (A > 100) a slightly different value is used
Q_HEAVY = 0.70


def effective_ga(q: float = Q_STANDARD) -> float:
    """Return effective axial coupling g_A^eff = q * g_A."""
    return q * G_A


def apply_gt_quenching(M_GT: float, q: float = Q_STANDARD) -> float:
    """
    Apply quenching to the GT matrix element.
        M_GT^eff = q * M_GT
    """
    return q * M_GT


def apply_ff_quenching(xi: float, zeta: float,
                       q_gt: float = Q_STANDARD,
                       q_ff: float = 1.0) -> tuple:
    """
    Apply quenching to first-forbidden matrix elements.

    GT-type (spin-dependent) part is quenched by q_gt.
    FF rank-0 part (xi) may be quenched differently.
    By default q_ff=1.0 (no quenching for FF operators).

    Returns
    -------
    tuple : (xi_eff, zeta_eff)
    """
    return q_ff * xi, q_gt * zeta


def quenching_from_mass(A: int) -> float:
    """
    Mass-number-dependent quenching factor.
    Simple empirical formula used for heavy nuclei.
    """
    if A <= 60:
        return Q_STANDARD
    elif A <= 100:
        return 0.72
    else:
        return Q_HEAVY
