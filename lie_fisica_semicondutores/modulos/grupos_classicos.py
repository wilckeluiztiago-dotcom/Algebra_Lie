"""
Módulo: grupos_classicos
Álgebras clássicas concretas usadas em física de partículas.

Referências no tratado:
- Cap. 1 §1.4.2  gl(n), sl(n)
- Cap. 1 §1.6.1  sl(2,C) – “o átomo da teoria”
- Cap. 1 §1.6.2  so(3,R) e momento angular
- Cap. 1 §1.6.3  Álgebra de Heisenberg h₃
- Cap. 2 §2.9    Zoológico de grupos clássicos
- Cap. 10 §10.2–10.3  Heisenberg e su(2) na QM
"""

from __future__ import annotations
import numpy as np
from .algebra_lie_basica import AlgebraDeLie


def geradores_pauli() -> List[np.ndarray]:
    """
    Matrizes de Pauli σ₁, σ₂, σ₃ (base hermitiana de su(2)).
    Relação: [σ_i/2, σ_j/2] = i ε_ijk (σ_k/2)
    """
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [sigma1, sigma2, sigma3]


def algebra_su2() -> AlgebraDeLie:
    """
    su(2) ≅ so(3) – álgebra do momento angular / isospin / spin.
    Constantes de estrutura: c^k_{ij} = ε_ijk  (totalmente antissimétricas)
    (tratado Tabela 1.1 e §1.6.1–1.6.2)
    """
    # ε_ijk
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = 1
    eps[1, 2, 0] = 1
    eps[2, 0, 1] = 1
    eps[0, 2, 1] = -1
    eps[2, 1, 0] = -1
    eps[1, 0, 2] = -1

    # Em física usa-se frequentemente geradores T_i = σ_i/2,
    # então [T_i, T_j] = i ε_ijk T_k  →  constantes complexas.
    # Aqui trabalhamos com a versão real so(3) (constantes reais).
    return AlgebraDeLie(
        dimensao=3,
        constantes_estrutura=eps,
        nome="su(2) ≅ so(3)",
        base_nomes=["J_x", "J_y", "J_z"],
    )


def algebra_so3() -> AlgebraDeLie:
    """Alias físico para a álgebra do momento angular clássico."""
    return algebra_su2()


def geradores_gell_mann() -> List[np.ndarray]:
    """
    Oito matrizes de Gell-Mann λ_a (a=1…8) – base hermitiana de su(3).
    Constantes de estrutura f_abc completamente antissimétricas
    (tratado Cap. 10 §10.6 e Apêndice B).
    """
    lambda1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    lambda2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    lambda3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    lambda4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    lambda5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    lambda6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    lambda7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    lambda8 = (1 / np.sqrt(3)) * np.array(
        [[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex
    )
    return [lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8]


def _constantes_estrutura_su3() -> np.ndarray:
    """
    Constantes de estrutura f_abc de su(3) (não nulas):
    f_123 = 1, f_147 = f_246 = f_257 = f_345 = 1/2,
    f_156 = f_367 = -1/2, f_458 = f_678 = √3/2
    (valores padrão da literatura, consistentes com o tratado).
    """
    f = np.zeros((8, 8, 8))
    # valores não nulos (índices 0-based)
    nao_nulos = {
        (0, 1, 2): 1.0,
        (0, 3, 6): 0.5,
        (0, 4, 5): -0.5,
        (1, 3, 5): 0.5,
        (1, 4, 6): 0.5,
        (2, 3, 4): 0.5,
        (2, 5, 6): -0.5,
        (3, 4, 7): np.sqrt(3) / 2,
        (5, 6, 7): np.sqrt(3) / 2,
    }
    for (i, j, k), val in nao_nulos.items():
        f[i, j, k] = val
        f[j, k, i] = val
        f[k, i, j] = val
        f[j, i, k] = -val
        f[i, k, j] = -val
        f[k, j, i] = -val
    return f


def algebra_su3() -> AlgebraDeLie:
    """
    su(3) – álgebra de cor (QCD) e de sabor (ócteto de mésons).
    dim = 8 → 8 glúons (tratado Cap. 10 eq. 10.25)
    """
    return AlgebraDeLie(
        dimensao=8,
        constantes_estrutura=_constantes_estrutura_su3(),
        nome="su(3)",
        base_nomes=[f"T{a+1}" for a in range(8)],
    )


def algebra_heisenberg() -> AlgebraDeLie:
    """
    Álgebra de Heisenberg h₃ (dim 3):
    base {x, y, z} com [x,y]=z e z central.
    (Cap. 1 §1.6.3 e Cap. 10 §10.2 – oscilador harmônico)
    """
    # índices: 0=x, 1=y, 2=z
    c = np.zeros((3, 3, 3))
    c[0, 1, 2] = 1.0   # [x,y] = z
    c[1, 0, 2] = -1.0  # antissimetria
    return AlgebraDeLie(
        dimensao=3,
        constantes_estrutura=c,
        nome="heisenberg h₃",
        base_nomes=["x", "y", "z"],
    )


# Necessário para type hints
from typing import List
