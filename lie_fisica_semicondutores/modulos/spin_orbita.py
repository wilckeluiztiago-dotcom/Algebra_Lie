"""
Módulo: spin_orbita
Acoplamentos Rashba e Dresselhaus como Hamiltonianos de su(2).

Referências no tratado de Luiz Tiago Wilcke (2026):
- Cap. 1 §1.6.1–1.6.2  su(2) e momento angular
- Cap. 10 §10.3        Momento angular e a álgebra su(2)
- Cap. 3 §3.7          Representações e Clebsch-Gordan (para composição)
"""

from __future__ import annotations
import numpy as np
from typing import Tuple, Dict, List, Optional
from .grupos_classicos import geradores_pauli, algebra_su2


class HamiltonianosSpinOrbita:
    """
    Hamiltoniano de spin-órbita linear em 2DEG:

        H = (ħ²k²)/(2m*)  +  α_R (σ_x k_y − σ_y k_x)  +  β_D (σ_x k_x − σ_y k_y)

    Os termos lineares são da forma Ω(k)·σ, com σ ∈ su(2).
    """

    def __init__(
        self,
        m_efetiva: float,
        alpha_rashba: float,
        beta_dresselhaus: float = 0.0,
        hbar: float = 1.0545718e-34,
    ):
        self.m_efetiva = m_efetiva
        self.alpha_R = alpha_rashba
        self.beta_D = beta_dresselhaus
        self.hbar = hbar
        self.sigma = geradores_pauli()  # σ_x, σ_y, σ_z

    def hamiltoniano_k(
        self, kx: float, ky: float
    ) -> np.ndarray:
        """
        Matriz 2×2 do Hamiltoniano em um ponto k do espaço recíproco.
        """
        k2 = kx**2 + ky**2
        H_cin = (self.hbar**2 * k2) / (2 * self.m_efetiva) * np.eye(2, dtype=complex)

        # Rashba: α (σ_x ky − σ_y kx)
        H_R = self.alpha_R * (
            self.sigma[0] * ky - self.sigma[1] * kx
        )

        # Dresselhaus linear: β (σ_x kx − σ_y ky)
        H_D = self.beta_D * (
            self.sigma[0] * kx - self.sigma[1] * ky
        )

        return H_cin + H_R + H_D

    def vetor_precessao(self, kx: float, ky: float) -> np.ndarray:
        """
        Vetor Ω(k) tal que H_SO = Ω·σ.
        Ω_x = α ky + β kx
        Ω_y = −α kx − β ky
        Ω_z = 0
        """
        Omega_x = self.alpha_R * ky + self.beta_D * kx
        Omega_y = -self.alpha_R * kx - self.beta_D * ky
        Omega_z = 0.0
        return np.array([Omega_x, Omega_y, Omega_z])

    def espectro(self, kx: float, ky: float) -> Tuple[float, float]:
        """
        Autovalores analíticos E± = ħ²k²/2m* ± |Ω(k)|.
        """
        H = self.hamiltoniano_k(kx, ky)
        evals = np.linalg.eigvalsh(H)
        return float(evals[0]), float(evals[1])

    def malha_espectro(
        self,
        k_max: float,
        n_pontos: int = 101,
    ) -> Dict:
        """
        Calcula a superfície E±(kx,ky) em uma malha.
        """
        k = np.linspace(-k_max, k_max, n_pontos)
        KX, KY = np.meshgrid(k, k)
        E_mais = np.zeros_like(KX)
        E_menos = np.zeros_like(KX)
        for i in range(n_pontos):
            for j in range(n_pontos):
                em, ep = self.espectro(KX[i, j], KY[i, j])
                E_menos[i, j] = em
                E_mais[i, j] = ep
        return {
            "kx": KX,
            "ky": KY,
            "E_menos": E_menos,
            "E_mais": E_mais,
            "referencia": "Wilcke 2026, Cap. 1 §1.6 e Cap. 10 §10.3",
        }

    def polarizacao_spin(self, kx: float, ky: float) -> np.ndarray:
        """
        Direção de polarização de spin do estado de menor energia
        (paralela a −Ω).
        """
        Omega = self.vetor_precessao(kx, ky)
        norma = np.linalg.norm(Omega)
        if norma < 1e-30:
            return np.array([0.0, 0.0, 1.0])
        return -Omega / norma


def espectro_rashba(
    k: float,
    m_efetiva: float,
    alpha_R: float,
    hbar: float = 1.0545718e-34,
) -> Tuple[float, float]:
    """
    Espectro analítico puro de Rashba (β=0):

        E±(k) = ħ²k²/(2m*) ± α_R k
    """
    E0 = (hbar**2 * k**2) / (2 * m_efetiva)
    return E0 - alpha_R * k, E0 + alpha_R * k


def espectro_rashba_dresselhaus(
    kx: float,
    ky: float,
    m_efetiva: float,
    alpha_R: float,
    beta_D: float,
    hbar: float = 1.0545718e-34,
) -> Tuple[float, float]:
    """
    Espectro completo Rashba+Dresselhaus.
    """
    modelo = HamiltonianosSpinOrbita(m_efetiva, alpha_R, beta_D, hbar)
    return modelo.espectro(kx, ky)


def angulo_precessao_spin(
    kx: float,
    ky: float,
    alpha_R: float,
    beta_D: float,
    comprimento: float,
) -> float:
    """
    Ângulo de precessão de spin ao longo de uma trajetória retilínea
    de comprimento L (efeito Datta-Das / transistor de spin).

    θ = (2 m* / ħ²) |Ω| L   (aproximação semiclassica)
    """
    Omega = np.array([
        alpha_R * ky + beta_D * kx,
        -alpha_R * kx - beta_D * ky,
        0.0,
    ])
    return 2 * np.linalg.norm(Omega) * comprimento   # em unidades ħ=1, m*=1
