"""
Half-life calculation for beta decay.

Combines phase space integrals and transition strengths to compute
partial and total half-lives, then compares with experimental values.

For allowed GT decay:
    t_{1/2} = D / (g_A^2 * f * B_GT)
    where D = 6289 s

For multiple transitions:
    1/t_{1/2} = sum_i 1/t_i

Reference: Suhonen (2007); Behrens & Bühring (1982).
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional

from .constants import D_CONST, G_A, KAPPA
from .phase_space import phase_space_integral
from .shape_factor import allowed_shape_factor, first_forbidden_shape_factor
from .quenching import effective_ga


@dataclass
class Transition:
    """Single beta decay transition."""
    label: str
    Q_MeV: float
    B_GT: float = 0.0
    B_F: float = 0.0
    xi: float = 0.0    # first-forbidden rank-0
    zeta: float = 0.0  # first-forbidden rank-1
    Z_daughter: int = 0
    g_A_eff: float = G_A


@dataclass
class HalfLifeResult:
    """Result of half-life calculation."""
    label: str
    f: float                # phase space integral
    B: float                # effective transition strength
    t_half_s: float         # partial half-life [seconds]
    t_half_ms: float        # partial half-life [milliseconds]
    log_ft: float           # log10(f * t_{1/2})
    exp_t_half_ms: Optional[float] = None  # experimental value [ms]

    def ratio(self) -> Optional[float]:
        """Calculated / Experimental ratio."""
        if self.exp_t_half_ms is None or self.exp_t_half_ms <= 0:
            return None
        return self.t_half_ms / self.exp_t_half_ms

    def __str__(self) -> str:
        exp_str = f"{self.exp_t_half_ms:.1f} ms" if self.exp_t_half_ms else "N/A"
        ratio_str = f"{self.ratio():.3f}" if self.ratio() else "N/A"
        return (f"{self.label}: t_1/2 = {self.t_half_ms:.2f} ms "
                f"(exp: {exp_str}, ratio: {ratio_str}), "
                f"log(ft) = {self.log_ft:.2f}")


def partial_half_life_gt(Z: int, Q_MeV: float, B_GT: float,
                          g_A_eff: float = G_A) -> float:
    """
    Partial half-life for a single allowed GT transition.

    t_{1/2} = D / (g_A_eff^2 * f * B_GT)  [seconds]
    """
    if B_GT <= 0 or Q_MeV <= 0:
        return float("inf")

    C = allowed_shape_factor(B_GT)
    f = phase_space_integral(Z, Q_MeV, C)
    if f <= 0:
        return float("inf")
    return D_CONST / (g_A_eff ** 2 * f)


def partial_half_life_ff(Z: int, Q_MeV: float, xi: float, zeta: float,
                          g_A_eff: float = G_A) -> float:
    """
    Partial half-life for a first-forbidden (non-unique) transition.

    t_{1/2} = D / (g_A_eff^2 * f_eff)
    where f_eff includes the shape factor C(W) = (xi - zeta*W)^2 + ...
    """
    if Q_MeV <= 0:
        return float("inf")

    C = first_forbidden_shape_factor(xi, zeta)
    f = phase_space_integral(Z, Q_MeV, C)
    if f <= 0:
        return float("inf")
    return D_CONST / (g_A_eff ** 2 * f)


def total_half_life(transitions: List[Transition]) -> float:
    """
    Total half-life from multiple transitions (decay rates add):
        1/t_{1/2} = sum_i 1/t_i
    """
    total_rate = 0.0
    for tr in transitions:
        if tr.B_GT > 0:
            t_i = partial_half_life_gt(tr.Z_daughter, tr.Q_MeV,
                                        tr.B_GT, tr.g_A_eff)
        elif abs(tr.xi) > 0 or abs(tr.zeta) > 0:
            t_i = partial_half_life_ff(tr.Z_daughter, tr.Q_MeV,
                                        tr.xi, tr.zeta, tr.g_A_eff)
        else:
            continue
        if math.isfinite(t_i) and t_i > 0:
            total_rate += 1.0 / t_i

    if total_rate <= 0:
        return float("inf")
    return 1.0 / total_rate


def compute_result(transition: Transition,
                   exp_t_half_ms: Optional[float] = None) -> HalfLifeResult:
    """
    Compute half-life result for a single transition with optional comparison.
    """
    if transition.B_GT > 0:
        C = allowed_shape_factor(transition.B_GT)
        B_eff = transition.B_GT
    else:
        C = first_forbidden_shape_factor(transition.xi, transition.zeta)
        B_eff = transition.xi ** 2 + transition.zeta ** 2

    f = phase_space_integral(transition.Z_daughter, transition.Q_MeV, C)
    t_half_s = D_CONST / (transition.g_A_eff ** 2 * f) if f > 0 else float("inf")
    t_half_ms = t_half_s * 1e3
    log_ft_val = math.log10(f * t_half_s) if (f > 0 and math.isfinite(t_half_s)) else float("inf")

    return HalfLifeResult(
        label=transition.label,
        f=f,
        B=B_eff,
        t_half_s=t_half_s,
        t_half_ms=t_half_ms,
        log_ft=log_ft_val,
        exp_t_half_ms=exp_t_half_ms,
    )


def compare_with_experiment(results: List[HalfLifeResult]) -> None:
    """Print a comparison table of calculated vs experimental half-lives."""
    print(f"{'Nucleus':<12} {'Calc [ms]':>12} {'Exp [ms]':>12} {'Ratio':>8} {'log(ft)':>9}")
    print("-" * 60)
    for r in results:
        exp_str = f"{r.exp_t_half_ms:.1f}" if r.exp_t_half_ms else "---"
        ratio_str = f"{r.ratio():.3f}" if r.ratio() else "---"
        print(f"{r.label:<12} {r.t_half_ms:>12.2f} {exp_str:>12} {ratio_str:>8} {r.log_ft:>9.2f}")
