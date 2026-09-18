"""
Modelo 1 — Oscilador harmônico quântico via álgebra de Heisenberg
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6.3 e Cap. 10 §10.2 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Tuple, Optional
from .grupos_classicos import algebra_heisenberg


class OsciladorHarmonico:
    """
    Oscilador harmônico unidimensional.

    Álgebra de Heisenberg: [x, p] = i ħ  (tratado Cap. 1 §1.6.3)
    Operadores de ladder: [a, a†] = 1
    H = ħ ω (a† a + 1/2)
    """

    def __init__(
        self,
        omega: float = 1.0,
        massa: float = 1.0,
        hbar: float = 1.0,
        n_corte: int = 30,
    ):
        self.omega = omega
        self.massa = massa
        self.hbar = hbar
        self.n_corte = n_corte
        self.dim = n_corte + 1
        self._construir_operadores()

    def _construir_operadores(self) -> None:
        n = np.arange(self.dim)
        # a |n⟩ = √n |n-1⟩
        self.a = np.zeros((self.dim, self.dim), dtype=complex)
        for i in range(1, self.dim):
            self.a[i - 1, i] = np.sqrt(i)
        self.a_dag = self.a.conj().T
        self.N = self.a_dag @ self.a
        self.H = self.hbar * self.omega * (self.N + 0.5 * np.eye(self.dim))

        # x e p a partir de a, a†
        fator_x = np.sqrt(self.hbar / (2 * self.massa * self.omega))
        fator_p = np.sqrt(self.hbar * self.massa * self.omega / 2)
        self.x = fator_x * (self.a + self.a_dag)
        self.p = -1j * fator_p * (self.a - self.a_dag)

    def espectro(self) -> np.ndarray:
        """E_n = ħ ω (n + 1/2)"""
        return np.array([self.hbar * self.omega * (n + 0.5) for n in range(self.dim)])

    def verificar_ccr(self) -> dict:
        """Verifica [x,p] = i ħ e [a, a†] = 1 na truncagem."""
        comut_xp = self.x @ self.p - self.p @ self.x
        comut_aa = self.a @ self.a_dag - self.a_dag @ self.a
        # na truncagem o último estado quebra levemente; olhamos o bloco principal
        bloco = slice(0, self.dim - 1)
        return {
            "max_|[x,p] - iħ|": float(np.max(np.abs(comut_xp[bloco, bloco] - 1j * self.hbar * np.eye(self.dim - 1)))),
            "max_|[a,a†] - 1|": float(np.max(np.abs(comut_aa[bloco, bloco] - np.eye(self.dim - 1)))),
            "referencia": "Wilcke 2026, Cap. 1 §1.6.3 e Cap. 10 §10.2",
        }

    def estado_numero(self, n: int) -> np.ndarray:
        psi = np.zeros(self.dim, dtype=complex)
        if 0 <= n <= self.n_corte:
            psi[n] = 1.0
        return psi

    def valor_esperado(self, operador: np.ndarray, estado: np.ndarray) -> complex:
        return estado.conj() @ operador @ estado


def espectro_oscilador(omega: float, n_max: int, hbar: float = 1.0) -> np.ndarray:
    return np.array([hbar * omega * (n + 0.5) for n in range(n_max + 1)])


def verificar_heisenberg_algebra() -> dict:
    """Confirma a estrutura abstrata h₃ do tratado."""
    h3 = algebra_heisenberg()
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])
    z = h3.colchete_lie(x, y)
    return {
        "algebra": str(h3),
        "[x,y]": z.real.tolist(),
        "relacao": "[x,y]=z (central) → [a,a†]=1 após quantização",
        "referencia": "Wilcke 2026, Cap. 1 §1.6.3",
    }
