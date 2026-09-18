"""
Modelo 5 — Efeito Zeeman (spin em campo magnético)
Autor: Luiz Tiago Wilcke
Fonte: Cap. 10 §10.3 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple
from .momento_angular_su2 import MomentoAngular
from .grupos_classicos import geradores_pauli


class HamiltonianoZeeman:
    """
    H_Z = − μ · B = (g μ_B / ħ) J · B

    Para spin 1/2: H = (g μ_B / 2) σ · B  (em unidades adequadas)
    """

    def __init__(
        self,
        B: np.ndarray,
        g: float = 2.0,
        mu_B: float = 1.0,
        j: float = 0.5,
        hbar: float = 1.0,
    ):
        self.B = np.asarray(B, dtype=float)
        self.g = g
        self.mu_B = mu_B
        self.j = j
        self.hbar = hbar
        self.ma = MomentoAngular(j, hbar)
        self.H = self._construir()

    def _construir(self) -> np.ndarray:
        # H = (g μ_B / ħ) (Jx Bx + Jy By + Jz Bz)
        fator = self.g * self.mu_B / self.hbar
        H = fator * (
            self.ma.Jx * self.B[0]
            + self.ma.Jy * self.B[1]
            + self.ma.Jz * self.B[2]
        )
        return H

    def espectro(self) -> np.ndarray:
        return np.linalg.eigvalsh(self.H)

    def estados_proprios(self) -> Tuple[np.ndarray, np.ndarray]:
        return np.linalg.eigh(self.H)


def espectro_zeeman(
    B_mod: float,
    g: float = 2.0,
    mu_B: float = 1.0,
    j: float = 0.5,
) -> np.ndarray:
    """
    Campo ao longo de z: E_m = g μ_B B m
    (para j=1/2: ± g μ_B B / 2)
    """
    m_vals = np.arange(-j, j + 0.1, 1.0)
    return g * mu_B * B_mod * m_vals


def hamiltoniano_zeeman_spin12(Bx: float, By: float, Bz: float, omega_L: float = 1.0) -> np.ndarray:
    """
    H = (ħ ω_L / 2) σ · n̂   com ω_L = g μ_B B / ħ
    """
    sig = geradores_pauli()
    B = np.array([Bx, By, Bz])
    norma = np.linalg.norm(B)
    if norma < 1e-15:
        return np.zeros((2, 2), dtype=complex)
    n = B / norma
    H = 0.5 * omega_L * norma * (n[0] * sig[0] + n[1] * sig[1] + n[2] * sig[2])
    return H
