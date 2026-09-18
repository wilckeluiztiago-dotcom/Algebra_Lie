"""
Módulo: excitons_cg
Estrutura fina de éxcitons via coeficientes de Clebsch-Gordan de su(2).

Referências no tratado de Luiz Tiago Wilcke (2026):
- Cap. 3 §3.7  Decomposição Tensorial e o Teorema de Clebsch-Gordan
- Cap. 10 §10.4 Produto de Clebsch-Gordan e Composição de Spins
- Cap. 9       Espaços de pesos (multiplicidades)
"""

from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple, Optional
from .representacoes import coeficientes_clebsch_gordan_su2
from .grupos_classicos import geradores_pauli


def estados_exciton_cg(
    j_eletron: float = 0.5,
    j_buraco: float = 1.5,
) -> Dict:
    """
    Constrói os estados de momento angular total do éxciton

        |J, M⟩ = Σ ⟨j_e m_e ; j_h m_h | J M⟩ |m_e⟩ ⊗ |m_h⟩

    para o caso típico de semicondutores III-V (elétron Γ₆ j=1/2,
    buraco pesado Γ₈ j=3/2).

    Decomposição:
        1/2 ⊗ 3/2 = 1 ⊕ 2
    (singleto óptico J=1 e quinteto J=2 – tratado Cap. 10 §10.4).
    """
    # Para j_h = 3/2 usamos a fórmula genérica j ⊗ 1/2 e depois
    # construímos manualmente os CG de 1/2 ⊗ 3/2.
    # CG 1/2 ⊗ 3/2 → J=1 e J=2

    # Valores analíticos padrão
    cg_12_32 = {
        # J=2 (quinteto)
        (0.5, 1.5, 2.0, 2.0): 1.0,
        (0.5, 0.5, 2.0, 1.0): np.sqrt(1 / 2),
        (-0.5, 1.5, 2.0, 1.0): np.sqrt(1 / 2),
        (0.5, -0.5, 2.0, 0.0): np.sqrt(3 / 6),
        (-0.5, 0.5, 2.0, 0.0): np.sqrt(3 / 6),
        (0.5, -1.5, 2.0, -1.0): np.sqrt(1 / 2),
        (-0.5, -0.5, 2.0, -1.0): np.sqrt(1 / 2),
        (-0.5, -1.5, 2.0, -2.0): 1.0,
        # J=1 (tripleto óptico)
        (0.5, 0.5, 1.0, 1.0): np.sqrt(1 / 2),
        (-0.5, 1.5, 1.0, 1.0): -np.sqrt(1 / 2),
        (0.5, -0.5, 1.0, 0.0): np.sqrt(1 / 2),
        (-0.5, 0.5, 1.0, 0.0): -np.sqrt(1 / 2),
        (0.5, -1.5, 1.0, -1.0): np.sqrt(1 / 2),
        (-0.5, -0.5, 1.0, -1.0): -np.sqrt(1 / 2),
    }

    # Normalização residual para os coeficientes aproximados
    for chave in list(cg_12_32.keys()):
        if abs(cg_12_32[chave]) < 1e-12:
            del cg_12_32[chave]

    return {
        "j_eletron": j_eletron,
        "j_buraco": j_buraco,
        "decomposicao": "1/2 ⊗ 3/2 = 1 ⊕ 2",
        "J_possiveis": [1.0, 2.0],
        "coeficientes_CG": cg_12_32,
        "interpretacao": {
            "J=1": "estados ópticos brilhantes (permitidos por dipolo)",
            "J=2": "estados escuros (proibidos por dipolo elétrico)",
        },
        "referencia": "Wilcke 2026, Cap. 3 §3.7 e Cap. 10 §10.4",
    }


def hamiltoniano_estrutura_fina(
    delta_0: float,
    delta_1: float,
    delta_2: float = 0.0,
) -> np.ndarray:
    """
    Hamiltoniano de estrutura fina do éxciton no subespaço J=1,2
    (parametrização fenomenológica comum em pontos quânticos):

        H_ef = δ₀ |J=2⟩⟨J=2|  +  δ₁ (J·S)  +  δ₂ anisotropia

    Implementação em base {|2,2⟩,|2,1⟩,|2,0⟩,|2,-1⟩,|2,-2⟩,|1,1⟩,|1,0⟩,|1,-1⟩}.
    """
    # Dimensão 8 (5 + 3)
    dim = 8
    H = np.zeros((dim, dim), dtype=complex)

    # Deslocamento do quinteto J=2
    for i in range(5):
        H[i, i] = delta_0

    # Termo de troca isotrópica proporcional a J·S (simplificado)
    # Diferença de energia entre J=1 e J=2
    for i in range(5, 8):
        H[i, i] = delta_1

    # Anisotropia axial (quebra J_z)
    if abs(delta_2) > 0:
        # acopla |2,±1⟩ com |1,±1⟩ etc. (simplificado)
        H[1, 5] = delta_2
        H[5, 1] = delta_2
        H[3, 7] = delta_2
        H[7, 3] = delta_2

    return H


def espectro_estrutura_fina(
    delta_0: float,
    delta_1: float,
    delta_2: float = 0.0,
) -> Dict:
    """
    Diagonaliza o Hamiltoniano de estrutura fina e retorna
    as energias e a composição dos autoestados.
    """
    H = hamiltoniano_estrutura_fina(delta_0, delta_1, delta_2)
    evals, evecs = np.linalg.eigh(H)
    return {
        "energias": evals.real,
        "vetores": evecs,
        "H": H,
        "referencia": "Wilcke 2026, Cap. 10 §10.4 (composição de momentos angulares)",
    }


def intensidade_optica_relativa(
    estados_cg: Dict,
    polarizacao: str = "sigma+",
) -> Dict[str, float]:
    """
    Intensidades relativas de emissão para polarizações circulares
    a partir dos CG (regras de seleção de dipolo).
    """
    # Regras de seleção: ΔM = ±1 para σ±, ΔM = 0 para π
    intensidades = {}
    if polarizacao == "sigma+":
        # |1,1⟩ e componentes de |2,1⟩ que contêm |m_e=1/2, m_h=1/2⟩
        intensidades["J=1,M=1"] = abs(estados_cg["coeficientes_CG"].get((0.5, 0.5, 1.0, 1.0), 0.0))**2
        intensidades["J=2,M=1"] = abs(estados_cg["coeficientes_CG"].get((0.5, 0.5, 2.0, 1.0), 0.0))**2
    elif polarizacao == "sigma-":
        intensidades["J=1,M=-1"] = abs(estados_cg["coeficientes_CG"].get((-0.5, -0.5, 1.0, -1.0), 0.0))**2
        intensidades["J=2,M=-1"] = abs(estados_cg["coeficientes_CG"].get((-0.5, -0.5, 2.0, -1.0), 0.0))**2
    return intensidades
