"""Unit tests for quenching module."""

import pytest
from beta_decay.quenching import (
    effective_ga, apply_gt_quenching, apply_ff_quenching,
    quenching_from_mass, Q_STANDARD, Q_HEAVY
)
from beta_decay.constants import G_A


class TestEffectiveGa:
    def test_default_quenching(self):
        g_eff = effective_ga()
        assert abs(g_eff - Q_STANDARD * G_A) < 1e-10

    def test_custom_quenching(self):
        q = 0.80
        g_eff = effective_ga(q)
        assert abs(g_eff - q * G_A) < 1e-10

    def test_no_quenching(self):
        g_eff = effective_ga(q=1.0)
        assert abs(g_eff - G_A) < 1e-10


class TestApplyGtQuenching:
    def test_reduces_matrix_element(self):
        M_GT = 2.0
        M_eff = apply_gt_quenching(M_GT, q=0.74)
        assert abs(M_eff - 0.74 * 2.0) < 1e-10

    def test_zero_input(self):
        assert apply_gt_quenching(0.0) == 0.0

    def test_sign_preserved(self):
        M_eff = apply_gt_quenching(-1.5, q=0.74)
        assert M_eff < 0


class TestApplyFfQuenching:
    def test_default_no_ff_quenching(self):
        xi, zeta = 1.0, 2.0
        xi_eff, zeta_eff = apply_ff_quenching(xi, zeta)
        assert abs(xi_eff - xi) < 1e-10        # q_ff=1 by default
        assert abs(zeta_eff - 0.74 * zeta) < 1e-10

    def test_both_quenched(self):
        xi_eff, zeta_eff = apply_ff_quenching(1.0, 1.0, q_gt=0.7, q_ff=0.9)
        assert abs(xi_eff - 0.9) < 1e-10
        assert abs(zeta_eff - 0.7) < 1e-10


class TestQuenchingFromMass:
    def test_light_nuclei(self):
        q = quenching_from_mass(40)
        assert q == Q_STANDARD

    def test_heavy_nuclei(self):
        q = quenching_from_mass(208)
        assert q == Q_HEAVY

    def test_intermediate(self):
        q = quenching_from_mass(80)
        assert 0.7 <= q <= 0.75

    def test_monotonic_trend(self):
        q_light = quenching_from_mass(40)
        q_heavy = quenching_from_mass(208)
        assert q_heavy <= q_light
