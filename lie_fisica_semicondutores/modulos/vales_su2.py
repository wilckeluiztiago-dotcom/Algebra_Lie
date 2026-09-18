"""
Módulo: vales_su2
Hamiltonianos efetivos de vale (grafeno, TMD) como álgebras su(2)/su(4).

Referências no tratado de Luiz Tiago Wilcke (2026):
- Cap. 1 §1.6.1  su(2) – o átomo da teoria
- Cap. 10 §10.3  su(2) na mecânica quântica
- Cap. 3         Representações e decomposição tensorial
"""

from __future__ import annotations
import numpy as np
from typing import Tuple, Dict, Optional
from .grupos_classicos import geradores_pauli


class HamiltonianosVale:
    """
    Hamiltoniano efetivo de dois vales (K e K′) em materiais 2D:

        H_τ = v_F (τ σ_x k_x + σ_y k_y) + Δ σ_z

    onde τ = ±1 rotula o vale. A estrutura é isomorfa a um su(2)
    de pseudospin (sub-redes A/B).
    """

    def __init__(
        self,
        velocidade_fermi: float,
        gap: float = 0.0,
    ):
        self.v_F = velocidade_fermi
        self.Delta = gap
        self.sigma = geradores_pauli()

    def hamiltoniano_vale(
        self,
        kx: float,
        ky: float,
        tau: int = +1,
        hbar: float = 1.0545718e-34,
    ) -> np.ndarray:
        """
        Matriz 2×2 no espaço de sub-rede para um vale fixo τ.
        Energia = ħ v_F |k|  (unidades SI).
        """
        H = hbar * self.v_F * (
            tau * self.sigma[0] * kx + self.sigma[1] * ky
        ) + self.Delta * self.sigma[2]
        return H

    def espectro(
        self,
        kx: float,
        ky: float,
        tau: int = +1,
    ) -> Tuple[float, float]:
        """
        E± = ± sqrt( (v_F k)² + Δ² )
        (independente de τ para o espectro de energia).
        """
        H = self.hamiltoniano_vale(kx, ky, tau)
        evals = np.linalg.eigvalsh(H)
        return float(evals[0]), float(evals[1])

    def hamiltoniano_quatro_componentes(
        self,
        kx: float,
        ky: float,
    ) -> np.ndarray:
        """
        Hamiltoniano 4×4 (vale ⊗ sub-rede) sem spin:

            H = diag( H_{τ=+1}, H_{τ=-1} )
        """
        H_K = self.hamiltoniano_vale(kx, ky, tau=+1)
        H_Kp = self.hamiltoniano_vale(kx, ky, tau=-1)
        H4 = np.zeros((4, 4), dtype=complex)
        H4[0:2, 0:2] = H_K
        H4[2:4, 2:4] = H_Kp
        return H4

    def hamiltoniano_spin_vale(
        self,
        kx: float,
        ky: float,
        lambda_SOC: float = 0.0,
    ) -> np.ndarray:
        """
        Modelo de TMD (MoS₂ etc.) com acoplamento spin-órbita de vale:

            H = v_F (τ σ_x k_x + σ_y k_y) + (Δ/2) σ_z + τ s_z λ_SOC σ_z

        (espaço: vale ⊗ spin ⊗ sub-rede → 8 componentes; aqui
        retornamos o bloco 4×4 para um spin fixo e depois o outro).
        """
        # Versão simplificada 4×4 (vale ⊗ sub-rede) com SOC efetivo
        H4 = self.hamiltoniano_quatro_componentes(kx, ky)
        # termo τ λ σ_z  (já incluído de forma diagonal por vale)
        H4[0:2, 0:2] += lambda_SOC * self.sigma[2]
        H4[2:4, 2:4] -= lambda_SOC * self.sigma[2]
        return H4


def espectro_grafeno_massivo(
    k: float,
    v_F: float,
    Delta: float,
    hbar: float = 1.0545718e-34,
) -> Tuple[float, float]:
    """
    Espectro analítico de grafeno com gap (ou BN hexagonal):

        E±(k) = ± sqrt( (ħ v_F k)² + Δ² )
    """
    E = np.sqrt((hbar * v_F * k) ** 2 + Delta ** 2)
    return -E, E


def densidade_estados_2d_dirac(
    energia: float,
    v_F: float,
    Delta: float = 0.0,
    g_s: int = 2,
    g_v: int = 2,
    hbar: float = 1.0545718e-34,
) -> float:
    """
    Densidade de estados por unidade de área para o cone de Dirac
    (com degenerescência de spin g_s e de vale g_v):

        D(E) = (g_s g_v |E|) / (2 π (ħ v_F)²)   para |E| > Δ
    """
    if abs(energia) < abs(Delta):
        return 0.0
    return (g_s * g_v * abs(energia)) / (2 * np.pi * (hbar * v_F) ** 2)


def operador_inversao_vale() -> np.ndarray:
    """
    Operador que troca os vales K ↔ K′ (matriz 4×4).
    Útil para analisar simetrias.
    """
    # troca os blocos 0:2 ↔ 2:4
    tau_x = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
    ], dtype=complex)
    return tau_x


def simetria_quiral_grafeno() -> Dict:
    """
    Operador de quiralidade γ₅ = σ_z ⊗ τ_0  (em grafeno sem gap).
    Anticomuta com o Hamiltoniano: {H, γ₅} = 0.
    """
    gamma5 = np.diag([1, -1, 1, -1]).astype(complex)
    return {
        "gamma5": gamma5,
        "propriedade": "{H, γ₅} = 0  (simetria quiral)",
        "referencia": "Wilcke 2026, Cap. 1 (estrutura algébrica de su(2))",
    }
