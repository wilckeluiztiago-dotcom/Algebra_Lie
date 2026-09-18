"""
Modelo 3 — Adição de momentos angulares (Clebsch-Gordan)
Autor: Luiz Tiago Wilcke
Fonte: Cap. 3 §3.7 e Cap. 10 §10.4 do tratado (eq. 10.13)
"""

from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple
from .representacoes import coeficientes_clebsch_gordan_su2
from .momento_angular_su2 import MomentoAngular


def adicionar_spins(j1: float, j2: float) -> Dict:
    """
    Decomposição V^{j1} ⊗ V^{j2} = ⊕_J V^J
    com J de |j1-j2| até j1+j2.
    """
    J_min = abs(j1 - j2)
    J_max = j1 + j2
    Js = []
    j = J_min
    while j <= J_max + 1e-9:
        Js.append(j)
        j += 1.0
    dim_produto = int(2 * j1 + 1) * int(2 * j2 + 1)
    dim_soma = sum(int(2 * J + 1) for J in Js)
    return {
        "j1": j1,
        "j2": j2,
        "J_possiveis": Js,
        "dim_produto": dim_produto,
        "dim_soma_irreps": dim_soma,
        "decomposicao": f"{j1} ⊗ {j2} = " + " ⊕ ".join(str(J) for J in Js),
        "referencia": "Wilcke 2026, Cap. 10 §10.4 eq. (10.13)",
    }


def base_acoplada(j1: float, j2: float) -> Dict:
    """
    Constrói os estados |J M⟩ na base produto via coeficientes CG.
    Implementação completa para o caso mais usado: 1/2 ⊗ 1/2.
    """
    if abs(j1 - 0.5) < 1e-9 and abs(j2 - 0.5) < 1e-9:
        cg = coeficientes_clebsch_gordan_su2(0.5, 0.5)
        # base produto: |↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩
        # estados acoplados:
        # |1,1⟩ = |↑↑⟩
        # |1,0⟩ = (|↑↓⟩ + |↓↑⟩)/√2
        # |1,-1⟩ = |↓↓⟩
        # |0,0⟩ = (|↑↓⟩ - |↓↑⟩)/√2
        estados = {
            (1.0, 1.0): np.array([1, 0, 0, 0], dtype=complex),
            (1.0, 0.0): np.array([0, 1 / np.sqrt(2), 1 / np.sqrt(2), 0], dtype=complex),
            (1.0, -1.0): np.array([0, 0, 0, 1], dtype=complex),
            (0.0, 0.0): np.array([0, 1 / np.sqrt(2), -1 / np.sqrt(2), 0], dtype=complex),
        }
        return {
            "base_produto": ["|↑↑⟩", "|↑↓⟩", "|↓↑⟩", "|↓↓⟩"],
            "estados_acoplados": estados,
            "cg_dict": cg,
            "referencia": "Wilcke 2026, Cap. 3 §3.7 e Cap. 10 §10.4",
        }
    raise NotImplementedError("Caso geral: use coeficientes_clebsch_gordan_su2 para j2=1/2")


def matriz_mudanca_base_12_12() -> np.ndarray:
    """
    Matriz unitária que leva a base produto na base |J M⟩
    (ordem: |1,1⟩, |1,0⟩, |1,-1⟩, |0,0⟩).
    """
    return np.array([
        [1, 0, 0, 0],
        [0, 1 / np.sqrt(2), 1 / np.sqrt(2), 0],
        [0, 0, 0, 1],
        [0, 1 / np.sqrt(2), -1 / np.sqrt(2), 0],
    ], dtype=complex).T  # colunas = novos vetores
