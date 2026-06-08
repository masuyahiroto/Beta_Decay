"""
Unit tests for half_life module.

Reference:
- Free neutron: t_{1/2} ≈ 614 s (experiment: 611.0 ± 1.0 s from PDG 2022)
- For GT dominant decays: log(ft) ~ 4-6 typical
"""

import pytest
import math
from beta_decay.half_life import (
    partial_half_life_gt, partial_half_life_ff,
    total_half_life, compute_result, compare_with_experiment,
    Transition, HalfLifeResult,
)
from beta_decay.constants import D_CONST, G_A


class TestPartialHalfLifeGt:
    def test_positive_for_valid_decay(self):
        t = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=1.0)
        assert t > 0 and math.isfinite(t)

    def test_zero_bgt_gives_inf(self):
        t = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=0.0)
        assert math.isinf(t)

    def test_zero_q_gives_inf(self):
        t = partial_half_life_gt(Z=20, Q_MeV=0.0, B_GT=1.0)
        assert math.isinf(t)

    def test_larger_bgt_shorter_lifetime(self):
        t1 = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=1.0)
        t2 = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=2.0)
        assert t2 < t1

    def test_larger_q_shorter_lifetime(self):
        t1 = partial_half_life_gt(Z=20, Q_MeV=3.0, B_GT=1.0)
        t2 = partial_half_life_gt(Z=20, Q_MeV=6.0, B_GT=1.0)
        assert t2 < t1

    def test_quenching_increases_lifetime(self):
        t_free = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=1.0, g_A_eff=G_A)
        t_quen = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=1.0, g_A_eff=0.74 * G_A)
        assert t_quen > t_free

    def test_neutron_free_decay_order_of_magnitude(self):
        """
        Free neutron: n -> p, Q=0.7824 MeV, Z=1.
        Using B_GT=1, g_A_eff=1 (pure GT, no Fermi) gives t ~ D/f ~ 3700 s.
        Actual neutron t_{1/2}=614 s includes Fermi + GT with correct couplings.
        Check that our simplified result is in the right ballpark (1000-8000 s).
        """
        t = partial_half_life_gt(Z=1, Q_MeV=0.7824, B_GT=1.0, g_A_eff=1.0)
        assert 1000 < t < 8000


class TestPartialHalfLifeFf:
    def test_positive_for_valid_decay(self):
        t = partial_half_life_ff(Z=82, Q_MeV=5.0, xi=1.0, zeta=0.5)
        assert t > 0 and math.isfinite(t)

    def test_zero_q_gives_inf(self):
        t = partial_half_life_ff(Z=82, Q_MeV=0.0, xi=1.0, zeta=0.5)
        assert math.isinf(t)

    def test_larger_matrix_element_shorter_lifetime(self):
        t1 = partial_half_life_ff(Z=82, Q_MeV=5.0, xi=1.0, zeta=0.0)
        t2 = partial_half_life_ff(Z=82, Q_MeV=5.0, xi=2.0, zeta=0.0)
        assert t2 < t1


class TestTotalHalfLife:
    def test_single_transition(self):
        tr = Transition(label="GT", Q_MeV=5.0, B_GT=1.0, Z_daughter=20)
        t_total = total_half_life([tr])
        t_partial = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=1.0)
        assert abs(t_total - t_partial) / t_partial < 1e-6

    def test_two_transitions_shorter_than_each(self):
        tr1 = Transition(label="GT1", Q_MeV=5.0, B_GT=1.0, Z_daughter=20)
        tr2 = Transition(label="GT2", Q_MeV=4.0, B_GT=0.5, Z_daughter=20)
        t_total = total_half_life([tr1, tr2])
        t1 = partial_half_life_gt(Z=20, Q_MeV=5.0, B_GT=1.0)
        t2 = partial_half_life_gt(Z=20, Q_MeV=4.0, B_GT=0.5)
        assert t_total < min(t1, t2)

    def test_empty_transitions(self):
        t = total_half_life([])
        assert math.isinf(t)

    def test_zero_bgt_transition_ignored(self):
        tr1 = Transition(label="GT", Q_MeV=5.0, B_GT=1.0, Z_daughter=20)
        tr2 = Transition(label="zero", Q_MeV=5.0, B_GT=0.0, Z_daughter=20)
        t_with = total_half_life([tr1, tr2])
        t_without = total_half_life([tr1])
        assert abs(t_with - t_without) < 1e-6


class TestComputeResult:
    def test_returns_result_object(self):
        tr = Transition(label="Pb208", Q_MeV=4.0, B_GT=0.5,
                        Z_daughter=82, g_A_eff=0.74 * G_A)
        result = compute_result(tr, exp_t_half_ms=100.0)
        assert isinstance(result, HalfLifeResult)
        assert result.t_half_ms > 0
        assert result.exp_t_half_ms == 100.0
        assert math.isfinite(result.log_ft)

    def test_ratio_computed(self):
        tr = Transition(label="test", Q_MeV=5.0, B_GT=1.0, Z_daughter=20)
        result = compute_result(tr, exp_t_half_ms=500.0)
        assert result.ratio() is not None
        assert result.ratio() > 0

    def test_ratio_none_without_exp(self):
        tr = Transition(label="test", Q_MeV=5.0, B_GT=1.0, Z_daughter=20)
        result = compute_result(tr, exp_t_half_ms=None)
        assert result.ratio() is None

    def test_t_half_units_consistent(self):
        tr = Transition(label="test", Q_MeV=5.0, B_GT=1.0, Z_daughter=20)
        result = compute_result(tr)
        assert abs(result.t_half_ms - result.t_half_s * 1e3) < 1e-6
