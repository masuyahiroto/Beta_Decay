"""
Phase space integrals for beta decay.

The integrated phase space factor is:

    f = integral_1^W0  C(W) * F(Z, W) * p * W * (W0 - W)^2  dW

where:
    W  = total electron energy / (m_e c^2)
    p  = electron momentum / (m_e c)  = sqrt(W^2 - 1)
    W0 = Q / (m_e c^2) + 1  (endpoint energy)
    F(Z, W) = Fermi function (Coulomb correction)
    C(W)    = shape factor

Reference: Behrens & Bühring (1982); Suhonen (2007) "From Nucleons to Nucleus".
"""

import numpy as np
from scipy import integrate
from .constants import ALPHA, ME_MEV
import math


def fermi_function(Z: int, W: float) -> float:
    """
    Fermi function F(Z, W) for beta-minus decay.

    Uses the point-nucleus approximation:
        eta = alpha * Z * W / p
        F(Z, W) = 2*pi*eta / (1 - exp(-2*pi*eta))

    Parameters
    ----------
    Z : int
        Proton number of the daughter nucleus.
    W : float
        Total electron energy in units of m_e c^2. Must be > 1.
    """
    if W <= 1.0:
        return 0.0
    p = math.sqrt(W * W - 1.0)
    eta = ALPHA * Z * W / p
    two_pi_eta = 2.0 * math.pi * eta
    # Guard against overflow
    if two_pi_eta > 500:
        return two_pi_eta
    return two_pi_eta / (1.0 - math.exp(-two_pi_eta))


def fermi_function_finite_size(Z: int, W: float, R_fm: float = 0.0) -> float:
    """
    Fermi function with finite nuclear size correction (L_0 factor).
    For heavy nuclei (Z~82), this correction is ~few percent.

    R_fm: nuclear radius in fm. If 0, uses R = 1.2 * A^{1/3} fm
          (A is estimated from Z for this function).

    Correction factor: L_0 ≈ 1 + (13/60)(alpha*Z)^2 - W*R/(10) * ...
    (simplified version)
    """
    F0 = fermi_function(Z, W)
    if R_fm <= 0.0:
        return F0
    # Finite size correction (leading order)
    R_natural = R_fm / 197.3269804  # convert fm -> 1/MeV -> natural units
    gamma = math.sqrt(1.0 - (ALPHA * Z) ** 2)
    L0 = 1.0 + (13.0 / 60.0) * (ALPHA * Z) ** 2 - W * R_natural * (41.0 / 15.0) * ALPHA * Z
    return F0 * max(L0, 0.5)


def phase_space_integral(Z: int, Q_MeV: float, shape_factor_func,
                         n_points: int = 2000,
                         use_finite_size: bool = False,
                         R_fm: float = 0.0) -> float:
    """
    Compute the phase space integral:
        f = integral_1^W0  C(W) * F(Z, W) * p * W * (W0-W)^2  dW

    Parameters
    ----------
    Z : int
        Proton number of daughter nucleus.
    Q_MeV : float
        Q-value of the decay in MeV.
    shape_factor_func : callable
        C(W) function of electron energy W.
    n_points : int
        Number of integration points.
    use_finite_size : bool
        Apply finite nuclear size correction.
    R_fm : float
        Nuclear radius in fm (used if use_finite_size=True).

    Returns
    -------
    float
        Dimensionless phase space integral f.
    """
    W0 = Q_MeV / ME_MEV + 1.0

    if W0 <= 1.0:
        return 0.0

    ff_func = fermi_function_finite_size if use_finite_size else \
              (lambda z, w: fermi_function(z, w))

    def integrand(W: float) -> float:
        if W <= 1.0 or W >= W0:
            return 0.0
        p = math.sqrt(W * W - 1.0)
        F = ff_func(Z, W)
        C = shape_factor_func(W)
        return C * F * p * W * (W0 - W) ** 2

    result, _ = integrate.quad(integrand, 1.0, W0, limit=200)
    return result


def log_ft(f: float, B: float, g_A: float = 1.2756) -> float:
    """
    Compute log(ft) value.
        log10(f * t_{1/2})

    Parameters
    ----------
    f    : float  - phase space integral
    B    : float  - transition strength B(GT) or B(F)
    g_A  : float  - axial coupling (effective)
    """
    from .constants import D_CONST
    if f <= 0 or B <= 0:
        return float("inf")
    t_half = D_CONST / (g_A ** 2 * f * B)
    return math.log10(f * t_half)
