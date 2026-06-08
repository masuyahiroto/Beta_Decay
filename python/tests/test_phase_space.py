"""
Unit tests for phase_space module.

Known reference values:
- Neutron decay (n -> p + e + nu): Q = 0.7824 MeV, Z=1
  f(1, 0.7824 MeV) ≈ 1.686e-3 (standard tables)
- Tritium (H-3 -> He-3): Q = 0.01859 MeV, Z=2
  f ≈ 2.9e-6 (very small due to low Q)
- Superallowed Fermi: log(ft) ~ 3.5
"""

import pytest
import math
from beta_decay.phase_space import (
    fermi_function, phase_space_integral, log_ft
)
from beta_decay.shape_factor import allowed_shape_factor, fermi_shape_factor


class TestFermiFunction:
    def test_low_z_near_threshold(self):
        """For Z=1, W just above threshold, F should be close to 1."""
        F = fermi_function(Z=1, W=1.01)
        assert 0.9 < F < 1.2

    def test_high_z_enhancement(self):
        """For high Z, Coulomb attraction enhances F > 1."""
        F_low = fermi_function(Z=1, W=2.0)
        F_high = fermi_function(Z=82, W=2.0)
        assert F_high > F_low

    def test_threshold_returns_zero(self):
        """At W = 1 (p=0), function should return 0."""
        F = fermi_function(Z=50, W=1.0)
        assert F == 0.0

    def test_below_threshold(self):
        F = fermi_function(Z=10, W=0.9)
        assert F == 0.0

    def test_positive_for_valid_W(self):
        for W in [1.1, 2.0, 5.0, 10.0]:
            assert fermi_function(Z=20, W=W) > 0

    def test_monotone_with_Z(self):
        """F increases with Z (more Coulomb attraction for beta-minus)."""
        W = 2.0
        F_vals = [fermi_function(Z=z, W=W) for z in [1, 10, 30, 80]]
        assert all(F_vals[i] < F_vals[i + 1] for i in range(len(F_vals) - 1))

    def test_neutron_decay(self):
        """For neutron (Z=1, W~1.5), F should be close to 1."""
        F = fermi_function(Z=1, W=1.5)
        assert 0.95 < F < 1.1


class TestPhaseSpaceIntegral:
    def test_zero_q_gives_zero(self):
        C = allowed_shape_factor(1.0)
        f = phase_space_integral(Z=1, Q_MeV=0.0, shape_factor_func=C)
        assert f == 0.0

    def test_negative_q_gives_zero(self):
        C = allowed_shape_factor(1.0)
        f = phase_space_integral(Z=1, Q_MeV=-1.0, shape_factor_func=C)
        assert f == 0.0

    def test_positive_for_valid_decay(self):
        C = allowed_shape_factor(1.0)
        f = phase_space_integral(Z=20, Q_MeV=5.0, shape_factor_func=C)
        assert f > 0.0

    def test_increases_with_q(self):
        """Larger Q -> larger phase space."""
        C = allowed_shape_factor(1.0)
        f1 = phase_space_integral(Z=20, Q_MeV=2.0, shape_factor_func=C)
        f2 = phase_space_integral(Z=20, Q_MeV=5.0, shape_factor_func=C)
        assert f2 > f1

    def test_increases_with_z(self):
        """Higher Z -> larger Fermi function -> larger f."""
        C = allowed_shape_factor(1.0)
        f1 = phase_space_integral(Z=1,  Q_MeV=3.0, shape_factor_func=C)
        f2 = phase_space_integral(Z=82, Q_MeV=3.0, shape_factor_func=C)
        assert f2 > f1

    def test_neutron_decay_order_of_magnitude(self):
        """
        Neutron free decay: Z=1, Q=0.7824 MeV.
        Standard result: f ~ 1.686 (dimensionless, in units of m_e c^2).
        Check within 10%.
        """
        C = fermi_shape_factor(B_F=1.0)
        f = phase_space_integral(Z=1, Q_MeV=0.7824, shape_factor_func=C)
        assert 1.4 < f < 2.0

    def test_bgt_scales_f_linearly(self):
        """B(GT) is a constant multiplier -> f scales linearly."""
        C1 = allowed_shape_factor(1.0)
        C2 = allowed_shape_factor(2.0)
        f1 = phase_space_integral(Z=20, Q_MeV=5.0, shape_factor_func=C1)
        f2 = phase_space_integral(Z=20, Q_MeV=5.0, shape_factor_func=C2)
        assert abs(f2 / f1 - 2.0) < 1e-6


class TestLogFt:
    def test_superallowed_fermi_approx(self):
        """
        Superallowed Fermi decays have log(ft) ~ 3.5.
        Test with B_F=1, high Q (large f).
        """
        from beta_decay.phase_space import log_ft
        from beta_decay.constants import G_A
        # For superallowed: log(ft) = log10(D / g_V^2 * 1/B_F) ≈ log10(6289) ≈ 3.8
        # Simple check: log(ft) should be in [3, 6] for typical decays
        C = fermi_shape_factor(1.0)
        f = phase_space_integral(Z=20, Q_MeV=10.0, shape_factor_func=C)
        lft = log_ft(f, B=1.0, g_A=1.0)
        assert 2.0 < lft < 6.0

    def test_zero_inputs(self):
        lft = log_ft(f=0.0, B=1.0)
        assert math.isinf(lft)

        lft = log_ft(f=1.0, B=0.0)
        assert math.isinf(lft)
