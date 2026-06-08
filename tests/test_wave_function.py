"""Unit tests for wave_function module."""

import pytest
from beta_decay.wave_function import (
    SingleParticleState, NuclearState, make_single_particle_transition
)


class TestSingleParticleState:
    def test_quantum_numbers(self):
        sp = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        assert sp.n == 0
        assert sp.l == 3
        assert sp.j == 3.5
        assert sp.tz == -0.5

    def test_repr_proton(self):
        sp = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        assert "f" in repr(sp)
        assert "p" in repr(sp)

    def test_repr_neutron(self):
        sp = SingleParticleState(n=0, l=3, j2=7, tz2=+1)
        assert "n" in repr(sp)

    def test_hashable(self):
        sp1 = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        sp2 = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        assert sp1 == sp2
        d = {sp1: 1.0}
        assert d[sp2] == 1.0


class TestNuclearState:
    def setup_method(self):
        self.p_state = SingleParticleState(n=0, l=3, j2=7, tz2=-1)  # 1f7/2 proton
        self.n_state = SingleParticleState(n=0, l=3, j2=7, tz2=+1)  # 1f7/2 neutron

    def test_set_get_obdm(self):
        state = NuclearState(J2=7, parity=-1, Z=20, N=20)
        state.set_obdm(self.p_state, self.n_state, 0.8)
        assert abs(state.get_obdm(self.p_state, self.n_state) - 0.8) < 1e-10

    def test_missing_obdm_returns_zero(self):
        state = NuclearState(J2=7, parity=-1, Z=20, N=20)
        assert state.get_obdm(self.p_state, self.n_state) == 0.0

    def test_mass_number(self):
        state = NuclearState(J2=7, parity=-1, Z=82, N=126)
        assert state.A == 208

    def test_spin(self):
        state = NuclearState(J2=7, parity=-1, Z=82, N=126)
        assert state.J == 3.5

    def test_len(self):
        state = NuclearState(J2=0, parity=+1, Z=8, N=8)
        state.set_obdm(self.p_state, self.n_state, 1.0)
        assert len(state) == 1


class TestMakeSingleParticleTransition:
    def test_creates_state_with_obdm(self):
        p = SingleParticleState(n=0, l=3, j2=7, tz2=-1)
        n = SingleParticleState(n=0, l=3, j2=7, tz2=+1)
        state = make_single_particle_transition(p, n, Z=20, N=21)
        assert state.get_obdm(p, n) == 1.0
        assert state.Z == 20
        assert state.N == 21
