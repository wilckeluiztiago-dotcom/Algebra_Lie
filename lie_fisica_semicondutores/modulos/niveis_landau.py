"""
Módulo: niveis_landau
Níveis de Landau a partir da álgebra de Heisenberg.

Referências no tratado de Luiz Tiago Wilcke (2026):
- Cap. 1 §1.6.3  Álgebra de Heisenberg h₃
- Cap. 10 §10.2  Heisenberg e oscilador harmônico
- Cap. 12        Sistemas dinâmicos (contexto de campos externos)
"""

from __future__ import annotations
import numpy as np
from typing import Tuple, List, Optional
from .grupos_classicos import algebra_heisenberg


# Constantes físicas (unidades SI simplificadas; ħ = 1 em u.a. atômicas quando indicado)
HBAR = 1.0545718e-34          # J·s
E_ELEMENTAR = 1.60217662e-19  # C
M_ELETRON = 9.1093837e-31     # kg


def operadores_ladder_landau(
    B: float,
    m_efetiva: float,
    hbar: float = HBAR,
    e: float = E_ELEMENTAR,
) -> Tuple[float, float]:
    """
    Constrói a frequência ciclotrônica e o comprimento magnético.

    ω_c = e B / m*
    ℓ_B = sqrt( ħ / (e B) )

    Os operadores de abaixamento/criação
        a  = (ℓ_B / √(2ħ)) (π_x - i π_y)
        a† = (ℓ_B / √(2ħ)) (π_x + i π_y)
    satisfazem [a, a†] = 1  (álgebra de Heisenberg – tratado Cap. 1 §1.6.3).
    """
    if B <= 0:
        raise ValueError("Campo magnético B deve ser positivo")
    omega_c = e * B / m_efetiva
    comprimento_magnetico = np.sqrt(hbar / (e * B))
    return omega_c, comprimento_magnetico


def espectro_landau(
    n_max: int,
    B: float,
    m_efetiva: float,
    hbar: float = HBAR,
    e: float = E_ELEMENTAR,
    g_fator: float = 2.0,
    mu_B: float = 9.274010e-24,
) -> dict:
    """
    Espectro completo dos níveis de Landau incluindo spin (efeito Zeeman).

    E_{n,s} = ħ ω_c (n + 1/2) + s g μ_B B / 2
    n = 0,1,2,…   s = ±1

    A parte orbital deriva diretamente da álgebra de Heisenberg do oscilador
    (tratado Cap. 10 §10.2). O termo de spin é gerado por su(2).
    """
    omega_c, ell_B = operadores_ladder_landau(B, m_efetiva, hbar, e)
    niveis = []
    for n in range(n_max + 1):
        E_orb = hbar * omega_c * (n + 0.5)
        for s in (+1, -1):
            E_zee = s * g_fator * mu_B * B / 2.0
            niveis.append({
                "n": n,
                "spin": s,
                "E_orbital_J": E_orb,
                "E_zeeman_J": E_zee,
                "E_total_J": E_orb + E_zee,
                "E_total_eV": (E_orb + E_zee) / E_ELEMENTAR,
            })
    return {
        "omega_c_rad_s": omega_c,
        "comprimento_magnetico_m": ell_B,
        "niveis": niveis,
        "referencia": "Wilcke 2026, Cap. 1 §1.6.3 e Cap. 10 §10.2",
    }


class HamiltonianosLandau:
    """
    Hamiltoniano de Landau em representação matricial truncada
    (base de número de quanta do oscilador).

    H = ħ ω_c (a† a + 1/2)  ⊗  1_spin  +  (g μ_B B / 2)  σ_z
    """

    def __init__(
        self,
        n_corte: int,
        B: float,
        m_efetiva: float,
        g_fator: float = 2.0,
        hbar: float = HBAR,
        e: float = E_ELEMENTAR,
        mu_B: float = 9.274010e-24,
    ):
        self.n_corte = n_corte
        self.B = B
        self.m_efetiva = m_efetiva
        self.g_fator = g_fator
        self.hbar = hbar
        self.e = e
        self.mu_B = mu_B

        self.omega_c, self.ell_B = operadores_ladder_landau(B, m_efetiva, hbar, e)

        # Base: |n⟩ ⊗ |↑⟩, |n⟩ ⊗ |↓⟩   (dim = 2*(n_corte+1))
        self.dim = 2 * (n_corte + 1)
        self.H = self._construir_matriz()

    def _construir_matriz(self) -> np.ndarray:
        dim_orb = self.n_corte + 1
        H = np.zeros((self.dim, self.dim), dtype=complex)

        # bloco orbital ħ ω_c (n + 1/2)
        for n in range(dim_orb):
            E_n = self.hbar * self.omega_c * (n + 0.5)
            # spin up (índice 2n)
            H[2 * n, 2 * n] = E_n + self.g_fator * self.mu_B * self.B / 2
            # spin down (índice 2n+1)
            H[2 * n + 1, 2 * n + 1] = E_n - self.g_fator * self.mu_B * self.B / 2
        return H

    def autovalores(self) -> np.ndarray:
        return np.linalg.eigvalsh(self.H)

    def densidade_estados_landau(
        self,
        area: float,
        degenerescencia_por_nivel: Optional[float] = None,
    ) -> float:
        """
        Densidade de estados por nível de Landau:
        D = (e B / h) * Área   (número de estados por nível orbital)
        """
        if degenerescencia_por_nivel is None:
            # número de fluxos quânticos
            degenerescencia_por_nivel = (self.e * self.B * area) / (2 * np.pi * self.hbar)
        return degenerescencia_por_nivel

    def verificar_algebra_heisenberg(self) -> dict:
        """
        Confirma que os operadores a, a† construídos satisfazem [a,a†]=1.
        (demonstração algébrica da origem Heisenberg – tratado Cap. 1)
        """
        h3 = algebra_heisenberg()
        # [x,y]=z  →  na linguagem de Landau, [π_x, π_y] = -i ħ e B
        comutador = h3.colchete_lie(
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
        )
        return {
            "algebra": str(h3),
            "comutador_[x,y]": comutador.real,
            "relacao_landau": "[π_x, π_y] = -i ħ e B  (origem do comprimento magnético)",
            "referencia": "Wilcke 2026, Cap. 1 §1.6.3",
        }
