"""
Nuclear wave function representation.

Stores one-body density matrix (OBDM) elements which are the key input
for computing transition matrix elements from shell model calculations.

OBDM: rho_{pn}^(J) = <J_f || c†_p c_n || J_i> / sqrt(2*J_i + 1)
"""

from dataclasses import dataclass, field
from typing import Dict, Tuple


@dataclass(frozen=True)
class SingleParticleState:
    """Single-particle state identified by quantum numbers (n, l, 2j, 2tz)."""
    n: int    # radial quantum number (0, 1, 2, ...)
    l: int    # orbital angular momentum
    j2: int   # 2 * total angular momentum j (integer)
    tz2: int  # 2 * isospin projection: -1 proton, +1 neutron

    @property
    def j(self) -> float:
        return self.j2 / 2

    @property
    def tz(self) -> float:
        return self.tz2 / 2

    def __repr__(self) -> str:
        l_labels = "spdfghik"
        tz_label = "n" if self.tz2 > 0 else "p"
        return f"{self.n + 1}{l_labels[self.l]}{self.j2}/2({tz_label})"


class NuclearState:
    """
    Nuclear state with quantum numbers and OBDM.

    The OBDM element rho[(proton_state, neutron_state)] corresponds to the
    reduced matrix element <f || c†_p c_n || i> / sqrt(2*Ji+1).
    """

    def __init__(self, J2: int, parity: int, Z: int, N: int):
        self.J2 = J2          # 2 * spin
        self.parity = parity  # +1 or -1
        self.Z = Z            # proton number
        self.N = N            # neutron number
        self._obdm: Dict[Tuple[SingleParticleState, SingleParticleState], float] = {}

    @property
    def J(self) -> float:
        return self.J2 / 2

    @property
    def A(self) -> int:
        return self.Z + self.N

    def set_obdm(self, p_state: SingleParticleState,
                 n_state: SingleParticleState, value: float) -> None:
        self._obdm[(p_state, n_state)] = value

    def get_obdm(self, p_state: SingleParticleState,
                 n_state: SingleParticleState) -> float:
        return self._obdm.get((p_state, n_state), 0.0)

    def items(self):
        return self._obdm.items()

    def __len__(self) -> int:
        return len(self._obdm)


def make_single_particle_transition(
    p_state: SingleParticleState,
    n_state: SingleParticleState,
    Z: int,
    N: int,
) -> NuclearState:
    """
    Convenience function: create a state with a single OBDM entry of 1.0.
    Useful for testing with pure single-particle transitions.
    """
    state = NuclearState(J2=n_state.j2, parity=(-1) ** n_state.l, Z=Z, N=N)
    state.set_obdm(p_state, n_state, 1.0)
    return state
