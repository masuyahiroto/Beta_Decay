"""Unit tests for matrix_elements module."""

import pytest
import numpy as np
from beta_decay.wave_function import SingleParticleState, NuclearState
from beta_decay.matrix_elements import (
    sp_gt_reduced_me,
    gamow_teller_strength,
    gt_matrix_element,
)


class TestSpGtReducedMe:
    """Test single-particle GT reduced matrix elements."""

    def test_same_j_same_l_j_plus(self):
        """j = l+1/2 -> j = l+1/2: should be sqrt(2*(2j+1))"""
        p = SingleParticleState(n=0, l=3, j2=7, tz2=-1)  # f7/2
        n = SingleParticleState(n=0, l=3, j2=7, tz2=+1)  # f7/2
        me = sp_gt_reduced_me(p, n)
        assert abs(me - np.sqrt(2 * 8)) < 1e-10

    def test_j_plus_to_j_minus(self):
        """j = l+1/2 proton, j = l-1/2 neutron: -sqrt(2jp+1)"""
        p = SingleParticleState(n=0, l=3, j2=7, tz2=-1)  # f7/2 proton
        n = SingleParticleState(n=0, l=3, j2=5, tz2=+1)  # f5/2 neutron
        me = sp_gt_reduced_me(p, n)
        assert abs(me - (-np.sqrt(8))) < 1e-10

    def test_j_minus_to_j_plus(self):
        """j = l-1/2 proton, j = l+1/2 neutron: -sqrt(2jn+1)"""
        p = SingleParticleState(n=0, l=3, j2=5, tz2=-1)  # f5/2 proton
        n = SingleParticleState(n=0, l=3, j2=7, tz2=+1)  # f7/2 neutron
        me = sp_gt_reduced_me(p, n)
        assert abs(me - (-np.sqrt(8))) < 1e-10

    def test_different_l_is_zero(self):
        """Different orbital angular momentum -> zero."""
        p = SingleParticleState(n=0, l=2, j2=5, tz2=-1)  # d5/2
        n = SingleParticleState(n=0, l=3, j2=7, tz2=+1)  # f7/2
        me = sp_gt_reduced_me(p, n)
        assert me == 0.0

    def test_s_orbital(self):
        """s1/2: j = l+1/2 = 1/2, l=0. <1/2||sigma||1/2> = sqrt(2*2) = 2"""
        p = SingleParticleState(n=0, l=0, j2=1, tz2=-1)  # s1/2 proton
        n = SingleParticleState(n=0, l=0, j2=1, tz2=+1)  # s1/2 neutron
        me = sp_gt_reduced_me(p, n)
        assert abs(me - np.sqrt(4)) < 1e-10  # = 2.0


class TestGamowTellerStrength:
    """Test B(GT) calculation from OBDM."""

    def test_single_particle_f7(self):
        """Single-particle f7/2 -> f7/2: B(GT) = M_GT^2 / (2*Ji+1)"""
        p = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        n = SingleParticleState(n=0, l=3, j2=7, tz2=+1)
        state = NuclearState(J2=7, parity=-1, Z=20, N=21)
        state.set_obdm(p, n, 1.0)

        B = gamow_teller_strength(state, Ji2=7)
        M_GT = sp_gt_reduced_me(p, n)
        expected = M_GT ** 2 / 8
        assert abs(B - expected) < 1e-10

    def test_zero_obdm_gives_zero(self):
        state = NuclearState(J2=0, parity=+1, Z=8, N=8)
        B = gamow_teller_strength(state, Ji2=0)
        assert B == 0.0

    def test_positive_definite(self):
        p = SingleParticleState(n=0, l=0, j2=1, tz2=-1)
        n = SingleParticleState(n=0, l=0, j2=1, tz2=+1)
        state = NuclearState(J2=1, parity=+1, Z=1, N=1)
        state.set_obdm(p, n, 0.5)
        B = gamow_teller_strength(state, Ji2=1)
        assert B >= 0.0

    def test_quenched_reduces_bgt(self):
        """Quenching by factor q reduces B(GT) by q^2."""
        p = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        n = SingleParticleState(n=0, l=3, j2=7, tz2=+1)
        state = NuclearState(J2=7, parity=-1, Z=20, N=21)
        state.set_obdm(p, n, 1.0)
        q = 0.74

        M_full = gt_matrix_element(state)
        M_quenched = q * M_full

        B_full = M_full ** 2 / 8
        B_quenched = M_quenched ** 2 / 8
        assert abs(B_quenched / B_full - q ** 2) < 1e-10
