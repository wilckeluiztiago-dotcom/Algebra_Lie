"""
Módulo: sistemas_raizes
Sistemas de raízes, matriz de Cartan e diagramas de Dynkin.

Referências no tratado:
- Capítulo 6: Sistemas de Raízes e Grupos de Weyl
- Capítulo 7: Matrizes de Cartan e Diagramas de Dynkin
- Teorema de classificação de Killing-Cartan (§7.2)
"""

from __future__ import annotations
import numpy as np
from typing import List, Dict, Tuple, Optional


class SistemaDeRaizes:
    """
    Sistema de raízes Φ ⊂ E (espaço euclidiano de posto ℓ).
    Axiomas (Cap. 6 §6.1):
    1. Φ finito, 0 ∉ Φ, span(Φ)=E
    2. α ∈ Φ ⇒ -α ∈ Φ
    3. α,β ∈ Φ ⇒ ⟨β,α∨⟩ ∈ Z  (número de Cartan)
    4. reflexão s_α(β) ∈ Φ
    """

    def __init__(
        self,
        raizes: np.ndarray,
        nome: str = "Φ",
    ):
        """
        raizes : array (N, ℓ) – lista de vetores-raiz
        """
        self.raizes = np.asarray(raizes, dtype=float)
        self.nome = nome
        self.posto = self.raizes.shape[1]
        self.numero_raizes = self.raizes.shape[0]

    def produto_interno(self, alpha: np.ndarray, beta: np.ndarray) -> float:
        return float(np.dot(alpha, beta))

    def numero_cartan(self, alpha: np.ndarray, beta: np.ndarray) -> int:
        """
        ⟨β, α∨⟩ = 2 (β·α) / (α·α)
        (definição Cap. 6 e 7 §7.1)
        """
        norma2 = self.produto_interno(alpha, alpha)
        if abs(norma2) < 1e-12:
            return 0
        return int(round(2 * self.produto_interno(beta, alpha) / norma2))

    def matriz_de_cartan(self, raizes_simples: np.ndarray) -> np.ndarray:
        """
        Matriz de Cartan A_{ij} = ⟨α_j, α_i∨⟩
        (Cap. 7 §7.1)
        """
        ell = raizes_simples.shape[0]
        A = np.zeros((ell, ell), dtype=int)
        for i in range(ell):
            for j in range(ell):
                A[i, j] = self.numero_cartan(raizes_simples[i], raizes_simples[j])
        return A


def matriz_de_cartan(tipo: str) -> np.ndarray:
    """
    Matrizes de Cartan das séries clássicas e de G₂ (exemplos do tratado).
    """
    tipo = tipo.upper()
    if tipo == "A1":
        return np.array([[2]])
    if tipo == "A2":
        return np.array([[2, -1], [-1, 2]])
    if tipo == "B2":
        return np.array([[2, -2], [-1, 2]])
    if tipo == "G2":
        # tratado Cap. 8 §8.2 – A_{21}=-3, A_{12}=-1
        return np.array([[2, -1], [-3, 2]])
    if tipo == "A3":
        return np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
    if tipo == "D4":
        return np.array(
            [
                [2, -1, 0, 0],
                [-1, 2, -1, -1],
                [0, -1, 2, 0],
                [0, -1, 0, 2],
            ]
        )
    raise ValueError(f"Tipo {tipo} ainda não implementado neste módulo.")


def diagrama_dynkin_simples(tipo: str) -> str:
    """
    Representação textual dos diagramas de Dynkin
    (Figuras 7.1 e 7.2 do tratado).
    """
    tipo = tipo.upper()
    diagramas = {
        "A1": "○",
        "A2": "○——○",
        "A3": "○——○——○",
        "B2": "○⇒⇒○  (raiz curta → longa)",
        "G2": "○≡≡≡○  (triplo, G₂)",
        "E6": "○——○——○——○——○\n        |\n        ○",
        "E7": "○——○——○——○——○——○\n        |\n        ○",
        "E8": "○——○——○——○——○——○——○\n        |\n        ○   (Fig. 8.3 – dim 248)",
        "F4": "○——○⇒⇒○——○",
    }
    return diagramas.get(tipo, f"(diagrama de {tipo} não listado)")


def sistema_raizes_a2() -> SistemaDeRaizes:
    """
    Sistema de raízes de A₂ = sl(3,C) (exemplo clássico Cap. 6 §6.4).
    6 raízes: ±α1, ±α2, ±(α1+α2)
    """
    alpha1 = np.array([1.0, 0.0])
    alpha2 = np.array([-0.5, np.sqrt(3) / 2])
    raizes = np.array(
        [
            alpha1,
            alpha2,
            alpha1 + alpha2,
            -alpha1,
            -alpha2,
            -(alpha1 + alpha2),
        ]
    )
    return SistemaDeRaizes(raizes, nome="A₂ = sl(3)")


def sistema_raizes_g2() -> SistemaDeRaizes:
    """
    Sistema de raízes de G₂ (Cap. 8 §8.2 e Fig. 8.1):
    6 curtas + 6 longas, razão √3 : 1.
    """
    # construção geométrica simplificada
    theta = np.pi / 6
    curtas = []
    longas = []
    for k in range(6):
        ang = k * np.pi / 3
        curtas.append([np.cos(ang), np.sin(ang)])
        longas.append([np.sqrt(3) * np.cos(ang + theta), np.sqrt(3) * np.sin(ang + theta)])
    raizes = np.array(curtas + longas)
    return SistemaDeRaizes(raizes, nome="G₂")
