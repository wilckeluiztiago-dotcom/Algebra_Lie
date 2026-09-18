"""
Modelo 10 — Parte angular do átomo de hidrogênio (so(3) / su(2))
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6.2 e Cap. 10 §10.3 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple
from scipy.special import sph_harm_y as sph_harm
from .momento_angular_su2 import MomentoAngular


class HarmonicosEsfericos:
    """
    Funções Y_lm — autovetores simultâneos de L² e L_z.

    L² Y_lm = ħ² l(l+1) Y_lm
    L_z Y_lm = ħ m Y_lm

    A álgebra é so(3) ≅ su(2) (tratado Cap. 1 §1.6.2).
    """

    def __init__(self, l_max: int = 3, hbar: float = 1.0):
        self.l_max = l_max
        self.hbar = hbar

    def Ylm(self, l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """
        Harmônico esférico Y_l^m (scipy.special.sph_harm_y(n, m, theta, phi)).
        """
        return sph_harm(l, m, theta, phi)

    def energia_angular_classica(self, l: int, r: float, massa: float = 1.0) -> float:
        """E_ang = ħ² l(l+1) / (2 m r²)"""
        return (self.hbar**2 * l * (l + 1)) / (2 * massa * r**2)

    def espectro_L2_Lz(self, l: int) -> Dict:
        ma = MomentoAngular(l, self.hbar)
        return {
            "l": l,
            "autovalor_L2": self.hbar**2 * l * (l + 1),
            "autovalores_Lz": ma.autovalores_jz(),
            "dim": 2 * l + 1,
            "referencia": "Wilcke 2026, Cap. 1 §1.6.2 e Cap. 10 §10.3",
        }

    def degenerescencia_nivel_n(self, n: int) -> int:
        """
        No hidrogênio, para o nível n a degenerescência (sem spin) é n²
        (soma de 2l+1 para l=0…n-1).
        """
        return n * n


def energia_angular(l: int, r: float, massa: float = 1.0, hbar: float = 1.0) -> float:
    return (hbar**2 * l * (l + 1)) / (2 * massa * r**2)


def tabela_degenerescencia(n_max: int = 4) -> Dict:
    """Tabela l, m e degenerescência até n = n_max."""
    dados = []
    for n in range(1, n_max + 1):
        for l in range(n):
            dados.append({
                "n": n,
                "l": l,
                "m_possiveis": list(range(-l, l + 1)),
                "degenerescencia_l": 2 * l + 1,
            })
    return {
        "niveis": dados,
        "degenerescencia_total_n": {n: n * n for n in range(1, n_max + 1)},
        "referencia": "Wilcke 2026, Cap. 1 §1.6.2 (so(3))",
    }
