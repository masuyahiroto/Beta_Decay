"""
Beta decay shape factor C(W).

The shape factor describes the W (electron energy) dependence of the
transition rate beyond the phase space factor.

Allowed (GT):
    C(W) = B_GT  (constant, no W dependence)

First-forbidden unique (Delta J = 2, no):
    C(W) = (W0 - W)^2 * |M|^2 + ...  (simplified)

First-forbidden non-unique (Delta J = 0, 1):
    C(W) = (xi - zeta * W)^2 + zeta^2 * (W^2 - 1) / 3  (rank decomposition)

Reference: Behrens & Bühring (1982); Kostensalo & Suhonen (2017).
"""

import numpy as np


def allowed_shape_factor(B_GT: float):
    """
    Shape factor for allowed Gamow-Teller transition.
    Returns a function C(W) = B_GT (constant).

    Parameters
    ----------
    B_GT : float
        Gamow-Teller transition strength.
    """
    def C(W: float) -> float:
        return B_GT

    return C


def first_forbidden_shape_factor(xi: float, zeta: float):
    """
    Shape factor for first-forbidden (Delta J = 0 or 1) transition.

    C(W) = (xi - zeta * W)^2 + zeta^2 * (W^2 - 1) / 3

    where:
        W   = total electron energy in units of m_e c^2
        xi  = rank-0 FF matrix element
        zeta= rank-1 FF matrix element

    Parameters
    ----------
    xi   : float
    zeta : float
    """
    def C(W: float) -> float:
        return (xi - zeta * W) ** 2 + zeta ** 2 * (W ** 2 - 1) / 3.0

    return C


def first_forbidden_unique_shape_factor(M: float):
    """
    Shape factor for first-forbidden unique (Delta J = 2) transition.

    C(W) = p^2 + (W0-W)^2  times |M|^2 (schematic)
    Here we use a simplified energy-independent form for integration.
    """
    def C(W: float) -> float:
        return M ** 2

    return C


def fermi_shape_factor(B_F: float):
    """Shape factor for allowed Fermi transition."""
    def C(W: float) -> float:
        return B_F

    return C
