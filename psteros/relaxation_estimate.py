"""How far an MLIP-relaxed structure is from its DFT minimum: a second-order estimate (no AiiDA).

On the MLIP route the model relaxes a structure and one DFT static at that
geometry gives its energy. The DFT forces ``F`` of that static (eV/A, on the
atoms that relax) and the model's Hessian ``H`` there (eV/A^2) give a
quasi-Newton step ``dx = H^-1 F`` towards the DFT minimum and the energy the
DFT relaxation would still lose,

``dE = -1/2 F^T H^-1 F``  (eV, <= 0),

exact for a quadratic energy surface. The modes of ``H`` that carry no
restoring force are left out: the ``zero_modes`` of smallest magnitude
(3 translations of a periodic cell whose every atom relaxes, 5 or 6 for a
molecule in a box) and every mode softer than ``min_curvature`` (a negative
curvature is a saddle point of the model, a very soft one a step the quadratic
model cannot be trusted with); the force on those modes is reported.

Applied to every structure of a study, ``E + dE`` give a second surface phase
diagram; the change of each gamma is the error bar of the route's geometry on
that termination (see :meth:`psteros.MlipRoute.gamma_error_bars`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence


@dataclass(frozen=True)
class RelaxationEstimate:
    """The second-order estimate of one structure (see :func:`relaxation_estimate`).

    ``energy_ev`` is ``dE`` (<= 0), ``max_force_ev_per_angstrom`` the largest
    DFT force on an atom that relaxes, ``max_step_angstrom`` the largest atom
    move of the step ``H^-1 F``, ``modes`` the modes used, ``negative_modes``
    and ``soft_modes`` those left out (negative curvature, or below
    ``min_curvature``), ``left_out_force_ev_per_angstrom`` the norm of the force
    on them. ``reliable`` is true when no mode was negative, the left-out force
    is below ``max_left_out_force`` and the step below ``max_step``.
    """

    energy_ev: float
    max_force_ev_per_angstrom: float
    max_step_angstrom: float
    modes: int
    negative_modes: int
    soft_modes: int
    left_out_force_ev_per_angstrom: float
    reliable: bool

    def summary(self) -> str:
        left_out = ""
        if self.negative_modes or self.soft_modes:
            left_out = (f", {self.negative_modes} negative and {self.soft_modes} soft modes left out "
                        f"(force {self.left_out_force_ev_per_angstrom:.3f} eV/A)")
        return (f"dE {1000 * self.energy_ev:8.2f} meV, largest DFT force {self.max_force_ev_per_angstrom:.3f} eV/A, "
                f"largest step {self.max_step_angstrom:.3f} A{left_out}{'' if self.reliable else ', NOT RELIABLE'}")

    def to_dict(self) -> dict[str, Any]:
        return {name: getattr(self, name) for name in self.__dataclass_fields__}


def relaxation_estimate(
    forces_ev_per_angstrom: Sequence[Sequence[float]],
    hessian_ev_per_angstrom2: Any,
    *,
    sites: Sequence[int],
    zero_modes: int = 0,
    min_curvature_ev_per_angstrom2: float = 0.05,
    max_step_angstrom: float = 0.1,
    max_left_out_force_ev_per_angstrom: float = 0.05,
) -> RelaxationEstimate:
    """``dE = -1/2 F^T H^-1 F`` over the modes of ``H`` with a restoring force.

    ``forces_ev_per_angstrom`` are the DFT forces on every atom of the
    structure (one row per site), ``hessian_ev_per_angstrom2`` the model's
    Hessian of ``sites`` (the atoms that relax, ``3 k + a`` order, e.g. from
    :func:`psteros.hessian_from_forces`). The ``zero_modes`` modes of smallest
    magnitude and every mode with curvature at most
    ``min_curvature_ev_per_angstrom2`` are left out. The estimate is
    ``reliable`` when no mode is negative, the force on the left-out modes is at
    most ``max_left_out_force_ev_per_angstrom`` and no atom moves more than
    ``max_step_angstrom``.
    """

    import numpy as np

    sites = [int(site) for site in sites]
    forces = np.asarray(forces_ev_per_angstrom, dtype=float)
    hessian = np.asarray(hessian_ev_per_angstrom2, dtype=float)
    if forces.ndim != 2 or forces.shape[1] != 3:
        raise ValueError("forces_ev_per_angstrom must have one row of 3 per site")
    if not sites or len(set(sites)) != len(sites) or any(site < 0 or site >= len(forces) for site in sites):
        raise ValueError(f"sites must be distinct site indices in [0, {len(forces)})")
    size = 3 * len(sites)
    if hessian.shape != (size, size):
        raise ValueError(f"hessian_ev_per_angstrom2 must have shape ({size}, {size}), got {hessian.shape}")
    if not np.allclose(hessian, hessian.T, atol=1e-8 * max(1.0, float(np.abs(hessian).max()))):
        raise ValueError("hessian_ev_per_angstrom2 must be symmetric")
    if not 0 <= zero_modes < size:
        raise ValueError(f"zero_modes must be in [0, {size})")
    if min_curvature_ev_per_angstrom2 <= 0 or max_step_angstrom <= 0 or max_left_out_force_ev_per_angstrom <= 0:
        raise ValueError("min_curvature, max_step and max_left_out_force must be positive")

    vector = forces[sites].reshape(size)
    curvatures, modes = np.linalg.eigh(hessian)
    zero = set(np.argsort(np.abs(curvatures))[:zero_modes].tolist())
    projections = modes.T @ vector
    used = [index for index in range(size) if index not in zero and curvatures[index] > min_curvature_ev_per_angstrom2]
    left_out = [index for index in range(size) if index not in zero and index not in used]
    negative = sum(1 for index in left_out if curvatures[index] < 0)
    step = modes[:, used] @ (projections[used] / curvatures[used]) if used else np.zeros(size)
    energy = -0.5 * float(np.sum(projections[used] ** 2 / curvatures[used])) if used else 0.0
    left_out_force = float(np.linalg.norm(projections[left_out])) if left_out else 0.0
    largest_step = float(np.linalg.norm(step.reshape(len(sites), 3), axis=1).max())
    return RelaxationEstimate(
        energy_ev=energy,
        max_force_ev_per_angstrom=float(np.linalg.norm(forces[sites], axis=1).max()),
        max_step_angstrom=largest_step,
        modes=len(used),
        negative_modes=negative,
        soft_modes=len(left_out) - negative,
        left_out_force_ev_per_angstrom=left_out_force,
        reliable=(negative == 0 and left_out_force <= max_left_out_force_ev_per_angstrom
                  and largest_step <= max_step_angstrom),
    )


@dataclass(frozen=True)
class GammaErrorBar:
    """How much the surface energy of one termination would change if DFT finished the relaxations.

    ``delta_gamma_ev_per_angstrom2`` is gamma(E + dE) - gamma(E) with every
    structure of the study (slabs, bulk, references) lowered by its estimate
    ``dE``, the largest change over the Delta mu range (it is the same at every
    Delta mu unless the references move differently); usually negative.
    ``estimate`` is that of the slab; ``reliable`` is true when the estimates of
    the slab and of every reference are reliable.
    """

    label: str
    delta_gamma_ev_per_angstrom2: float
    estimate: RelaxationEstimate
    reliable: bool

    def summary(self) -> str:
        return (f"{self.label}: d gamma {1000 * self.delta_gamma_ev_per_angstrom2:+.3f} meV/A2 "
                f"({self.estimate.summary()}){'' if self.reliable else ' [not reliable]'}")


def gamma_error_bars(base: Any, lowered: Any, estimates: Any, slabs: Sequence[str]) -> dict[str, GammaErrorBar]:
    """The :class:`GammaErrorBar` of each slab from two results of one study: ``base`` with the energies ``E`` and
    ``lowered`` with ``E + dE`` (their ``diagram.curves``), and the estimates ``{label: RelaxationEstimate}``."""

    references = [label for label in estimates if label not in set(slabs)]
    references_reliable = all(estimates[label].reliable for label in references)
    bars = {}
    for label in slabs:
        shifts = [after.gamma_ev_per_angstrom2 - before.gamma_ev_per_angstrom2
                  for before, after in zip(base.diagram.curves[label], lowered.diagram.curves[label])]
        largest = max(shifts, key=abs)
        bars[label] = GammaErrorBar(label, float(largest), estimates[label],
                                    estimates[label].reliable and references_reliable)
    return bars
