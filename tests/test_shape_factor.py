"""Unit tests for shape_factor module."""

import pytest
from beta_decay.shape_factor import (
    allowed_shape_factor,
    first_forbidden_shape_factor,
    first_forbidden_unique_shape_factor,
    fermi_shape_factor,
)


class TestAllowedShapeFactor:
    def test_constant_value(self):
        C = allowed_shape_factor(B_GT=2.5)
        assert C(1.5) == 2.5
        assert C(2.0) == 2.5
        assert C(5.0) == 2.5

    def test_energy_independent(self):
        C = allowed_shape_factor(B_GT=1.0)
        values = [C(W) for W in [1.1, 2.0, 3.5, 10.0]]
        assert all(v == 1.0 for v in values)

    def test_zero_bgt(self):
        C = allowed_shape_factor(B_GT=0.0)
        assert C(2.0) == 0.0


class TestFirstForbiddenShapeFactor:
    def test_positive_definite_for_real_inputs(self):
        C = first_forbidden_shape_factor(xi=1.0, zeta=0.5)
        for W in [1.1, 2.0, 5.0, 10.0]:
            assert C(W) >= 0.0

    def test_energy_dependence(self):
        """C(W) should vary with W when zeta != 0."""
        C = first_forbidden_shape_factor(xi=1.0, zeta=1.0)
        assert C(1.5) != C(3.0)

    def test_xi_only(self):
        """With zeta=0: C(W) = xi^2 (constant)."""
        xi = 1.5
        C = first_forbidden_shape_factor(xi=xi, zeta=0.0)
        assert abs(C(2.0) - xi ** 2) < 1e-10
        assert abs(C(5.0) - xi ** 2) < 1e-10

    def test_zeta_only(self):
        """With xi=0: C(W) = zeta^2*W^2 + zeta^2*(W^2-1)/3"""
        zeta = 1.0
        C = first_forbidden_shape_factor(xi=0.0, zeta=zeta)
        W = 2.0
        expected = (zeta * W) ** 2 + zeta ** 2 * (W ** 2 - 1) / 3.0
        assert abs(C(W) - expected) < 1e-10

    def test_formula_at_specific_point(self):
        xi, zeta, W = 1.0, 0.5, 2.0
        C = first_forbidden_shape_factor(xi=xi, zeta=zeta)
        expected = (xi - zeta * W) ** 2 + zeta ** 2 * (W ** 2 - 1) / 3.0
        assert abs(C(W) - expected) < 1e-10


class TestFermiShapeFactor:
    def test_constant(self):
        C = fermi_shape_factor(B_F=1.0)
        assert C(1.5) == 1.0
        assert C(100.0) == 1.0
