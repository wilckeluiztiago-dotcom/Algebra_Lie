"""
Modelo 6 — Sistema de dois níveis e precessão de Rabi
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6 e Cap. 10 (estrutura su(2))
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple, Callable
from scipy.linalg import expm
from .grupos_classicos import geradores_pauli


class HamiltonianoRabi:
    """
    Hamiltoniano de Rabi (dois níveis acoplados a um campo clássico):

    H = (ħ ω₀ / 2) σ_z + (ħ Ω / 2) (σ_x cos(ω t) + σ_y sin(ω t))

    Na aproximação de onda rotante (RWA) e no frame rotativo:
    H_RWA = (ħ Δ / 2) σ_z + (ħ Ω / 2) σ_x
    onde Δ = ω₀ − ω (dessintonia).
    """

    def __init__(
        self,
        omega_0: float,
        Omega: float,
        omega: float = None,
        hbar: float = 1.0,
    ):
        self.omega_0 = omega_0
        self.Omega = Omega
        self.omega = omega if omega is not None else omega_0
        self.hbar = hbar
        self.Delta = omega_0 - self.omega
        self.sigma = geradores_pauli()

    def H_lab(self, t: float) -> np.ndarray:
        """Hamiltoniano no laboratório."""
        H0 = 0.5 * self.hbar * self.omega_0 * self.sigma[2]
        H_int = 0.5 * self.hbar * self.Omega * (
            self.sigma[0] * np.cos(self.omega * t)
            + self.sigma[1] * np.sin(self.omega * t)
        )
        return H0 + H_int

    def H_rwa(self) -> np.ndarray:
        """Hamiltoniano efetivo na aproximação de onda rotante (frame rotativo)."""
        return 0.5 * self.hbar * (
            self.Delta * self.sigma[2] + self.Omega * self.sigma[0]
        )

    def frequencia_rabi_generalizada(self) -> float:
        """Ω_R = √(Ω² + Δ²)"""
        return np.sqrt(self.Omega**2 + self.Delta**2)

    def evolucao_rwa(self, t: float, psi0: np.ndarray) -> np.ndarray:
        """Evolução exata sob H_RWA."""
        U = expm(-1j * self.H_rwa() * t / self.hbar)
        return U @ psi0

    def probabilidade_excitacao(self, t: float, psi0: np.ndarray = None) -> float:
        """Probabilidade de encontrar o sistema no estado |e⟩ = |↑⟩."""
        if psi0 is None:
            psi0 = np.array([0.0, 1.0], dtype=complex)  # começa no |g⟩
        psi_t = self.evolucao_rwa(t, psi0)
        return float(np.abs(psi_t[0])**2)


def frequencia_rabi(Omega: float, Delta: float = 0.0) -> float:
    return np.sqrt(Omega**2 + Delta**2)


def oscilacoes_rabi(
    Omega: float,
    Delta: float,
    t_max: float,
    n_pontos: int = 200,
) -> Dict:
    """
    Curva P_e(t) para estado inicial |g⟩.
    P_e(t) = (Ω² / Ω_R²) sin²(Ω_R t / 2)
    """
    Omega_R = frequencia_rabi(Omega, Delta)
    t = np.linspace(0, t_max, n_pontos)
    if Omega_R < 1e-15:
        P = np.zeros_like(t)
    else:
        P = (Omega**2 / Omega_R**2) * np.sin(Omega_R * t / 2)**2
    return {
        "t": t,
        "P_excitacao": P,
        "Omega_R": Omega_R,
        "referencia": "Wilcke 2026, Cap. 1 e Cap. 10 (su(2))",
    }
