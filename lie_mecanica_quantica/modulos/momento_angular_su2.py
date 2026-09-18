"""
Modelo 2 — Momento angular e álgebra su(2)
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6.1–1.6.2 e Cap. 10 §10.3 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Tuple, List
from .grupos_classicos import algebra_su2, geradores_pauli


class MomentoAngular:
    """
    Representação irredutível de spin j de su(2).

    [J_i, J_j] = i ħ ε_ijk J_k
    J² |j m⟩ = ħ² j(j+1) |j m⟩
    J_z |j m⟩ = ħ m |j m⟩
    """

    def __init__(self, j: float, hbar: float = 1.0):
        if (2 * j) != int(2 * j) or j < 0:
            raise ValueError("j deve ser inteiro ou semi-inteiro não-negativo")
        self.j = j
        self.hbar = hbar
        self.dim = int(2 * j + 1)
        self.m_valores = np.array([j - k for k in range(self.dim)])
        self._construir_matrizes()

    def _construir_matrizes(self) -> None:
        j, h = self.j, self.hbar
        dim = self.dim
        self.Jz = h * np.diag(self.m_valores).astype(complex)

        # J± = Jx ± i Jy
        Jp = np.zeros((dim, dim), dtype=complex)
        Jm = np.zeros((dim, dim), dtype=complex)
        for i, m in enumerate(self.m_valores):
            # J+ |j m⟩ = ħ √(j(j+1)-m(m+1)) |j m+1⟩
            if m < j:
                idx_up = i - 1  # m aumenta → índice diminui na ordenação
                # reordenar: m = j, j-1, ..., -j  → índice 0,1,...,2j
                k = int(j - m)  # índice de m
                k_up = int(j - (m + 1))
                coef = h * np.sqrt(j * (j + 1) - m * (m + 1))
                Jp[k_up, k] = coef
                Jm[k, k_up] = coef  # J- = (J+)†

        self.Jp = Jp
        self.Jm = Jm
        self.Jx = 0.5 * (Jp + Jm)
        self.Jy = -0.5j * (Jp - Jm)
        self.J2 = self.Jx @ self.Jx + self.Jy @ self.Jy + self.Jz @ self.Jz

    def autovalores_j2(self) -> np.ndarray:
        return np.full(self.dim, self.hbar**2 * self.j * (self.j + 1))

    def autovalores_jz(self) -> np.ndarray:
        return self.hbar * self.m_valores

    def verificar_algebra(self) -> dict:
        """Verifica [Jx,Jy] = i ħ Jz etc."""
        h = self.hbar
        comut_xy = self.Jx @ self.Jy - self.Jy @ self.Jx
        esperado = 1j * h * self.Jz
        erro = np.max(np.abs(comut_xy - esperado))
        return {
            "max_|[Jx,Jy] - iħ Jz|": float(erro),
            "J2 autovalor teórico": self.hbar**2 * self.j * (self.j + 1),
            "referencia": "Wilcke 2026, Cap. 1 §1.6.1 e Cap. 10 §10.3",
        }

    def estado_jm(self, m: float) -> np.ndarray:
        psi = np.zeros(self.dim, dtype=complex)
        idx = list(self.m_valores).index(m)
        psi[idx] = 1.0
        return psi


def autovalores_j2_jz(j: float, hbar: float = 1.0) -> Tuple[float, np.ndarray]:
    """Retorna (j(j+1) ħ², lista de m ħ)."""
    m = np.arange(-j, j + 0.1, 1.0)
    return hbar**2 * j * (j + 1), hbar * m


def geradores_spin_1_2() -> List[np.ndarray]:
    """σ_i / 2  (geradores de su(2) na representação fundamental)."""
    return [s / 2 for s in geradores_pauli()]
