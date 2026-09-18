"""
Módulo: fisica_particulas
Aplicações diretas à física de partículas a partir do tratado.

Referências centrais:
- Capítulo 10: Mecânica Quântica e Física de Partículas
  §10.5 Teoria de Gauge e Equações de Yang-Mills
  §10.6 O Modelo Padrão: SU(3)×SU(2)×U(1)
  §10.6.1 Decomposição Tensorial: Mésons e Bárions
  §10.7 Teorias de Grande Unificação
- Capítulo 8: Álgebras Excepcionais (E₆, E₈ e GUTs)
"""

from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple
from .grupos_classicos import algebra_su2, algebra_su3, algebra_heisenberg
from .representacoes import dimensao_irrep_su3, coeficientes_clebsch_gordan_su2
from .algebra_lie_basica import AlgebraDeLie


def modelo_padrao() -> Dict:
    """
    Conteúdo de geradores do Modelo Padrão
    (tratado Cap. 10 §10.6, eqs 10.25–10.27):

        dim su(3)_c = 8  →  8 glúons
        dim su(2)_L = 3  →  W±, Z⁰
        dim u(1)_Y  = 1  →  fóton (após mistura de Weinberg)

    Total de bósons de gauge: 12.
    """
    su3 = algebra_su3()
    su2 = algebra_su2()
    return {
        "grupo": "SU(3)_c × SU(2)_L × U(1)_Y",
        "dim_su3": su3.dimensao,
        "dim_su2": su2.dimensao,
        "dim_u1": 1,
        "total_bosons_gauge": su3.dimensao + su2.dimensao + 1,
        "gluons": 8,
        "bosons_fracos": 3,  # W+, W-, Z
        "foton": 1,
        "referencia": "Wilcke 2026, Cap. 10 §10.6 eqs. (10.25)–(10.27)",
    }


def decomposicao_mesons() -> Dict:
    """
    Identificação de mésons via produto tensorial de sabor SU(3)
    (tratado eq. 10.28 e Fig. 10.2):

        3 ⊗ 3̄ = 1 ⊕ 8

    O octeto contém π, K, η; o singleto é o η'.
    """
    dim3 = dimensao_irrep_su3(1, 0)
    dim3bar = dimensao_irrep_su3(0, 1)
    dim1 = dimensao_irrep_su3(0, 0)
    dim8 = dimensao_irrep_su3(1, 1)
    assert dim3 * dim3bar == dim1 + dim8
    return {
        "produto": "3 ⊗ 3̄",
        "decomposicao": "1 ⊕ 8",
        "dimensoes": {"3": dim3, "3̄": dim3bar, "1": dim1, "8": dim8},
        "particulas_octeto": [
            "π⁺", "π⁰", "π⁻",
            "K⁺", "K⁰", "K̄⁰", "K⁻",
            "η",
        ],
        "referencia": "Wilcke 2026, Cap. 10 §10.6.1 eq. (10.28) e Fig. 10.2",
    }


def decomposicao_barions() -> Dict:
    """
    Bárions (tratado eq. 10.29):

        3 ⊗ 3 ⊗ 3 = 1 ⊕ 8 ⊕ 8 ⊕ 10

    O decupleto 10 contém o Ω⁻ (S=-3, spin 3/2), previsto e descoberto em 1964.
    """
    dim3 = dimensao_irrep_su3(1, 0)
    dim1 = dimensao_irrep_su3(0, 0)
    dim8 = dimensao_irrep_su3(1, 1)
    dim10 = dimensao_irrep_su3(3, 0)
    assert dim3 ** 3 == dim1 + 2 * dim8 + dim10
    return {
        "produto": "3 ⊗ 3 ⊗ 3",
        "decomposicao": "1 ⊕ 8 ⊕ 8 ⊕ 10",
        "dimensoes": {
            "3": dim3,
            "1": dim1,
            "8": dim8,
            "10": dim10,
        },
        "destaque": "Ω⁻ previsto no decupleto 10 (estranheza −3, spin 3/2)",
        "referencia": "Wilcke 2026, Cap. 10 §10.6.1 eq. (10.29)",
    }


def grupo_gut_su5() -> Dict:
    """
    SU(5) de Georgi-Glashow (tratado Cap. 10 §10.7 eq. 10.30):

        dim su(5) = 24 = 8 (glúons) + 3 (W,Z) + 1 (γ) + 12 (X,Y)

    Os 12 bósons X,Y medeiam o decaimento do próton p → e⁺ π⁰.
    """
    dim_su5 = 24  # n²-1 para n=5
    return {
        "grupo": "SU(5)",
        "dimensao": dim_su5,
        "conteudo": {
            "gluons": 8,
            "W_Z": 3,
            "foton": 1,
            "bosons_X_Y": 12,
        },
        "decaimento_proton": "p → e⁺ + π⁰  (via troca de X)",
        "tempo_vida_previsto": "~10³⁰ anos (excluído por Super-Kamiokande)",
        "referencia": "Wilcke 2026, Cap. 10 §10.7 eqs. (10.30)–(10.31)",
    }


def cadeia_excepcional_e8() -> Dict:
    """
    Cadeia de subgrupos dominante (tratado Cap. 8 §8.6 e Cap. 10 eq. 10.32):

        E₈ ⊃ E₇ ⊃ E₆ ⊃ SO(10) ⊃ SU(5) ⊃ SU(3)×SU(2)×U(1)

    Dados estruturais das excepcionais (Tabela 8.1 do tratado).
    """
    return {
        "cadeia": "E₈ ⊃ E₇ ⊃ E₆ ⊃ SO(10) ⊃ SU(5) ⊃ SU(3)×SU(2)×U(1)",
        "dimensoes": {
            "E₈": 248,
            "E₇": 133,
            "E₆": 78,
            "SO(10)": 45,
            "SU(5)": 24,
            "SM": 12,
        },
        "rep_fundamental_E6": {
            "27": "16 ⊕ 10 ⊕ 1  (sob SO(10))",
            "16": "spinor de SO(10) = uma geração completa de férmions (com ν_R)",
        },
        "referencia": "Wilcke 2026, Cap. 8 §8.5–8.6 e Cap. 10 §10.7",
    }


def tensor_yang_mills(
    potencial_A: np.ndarray,
    derivadas_parciais: np.ndarray,
    constante_acoplamento: float,
    constantes_estrutura: np.ndarray,
) -> np.ndarray:
    """
    Tensor de campo (curvatura) de Yang-Mills (tratado eq. 10.20):

        F_μν = ∂_μ A_ν − ∂_ν A_μ − i g [A_μ, A_ν]

    Implementação simbólica/matricial para ilustração.
    potencial_A : array (4, dim_g)  – componentes A_μ^a
    derivadas_parciais : array (4, 4, dim_g)  – ∂_μ A_ν^a
    """
    dim_espaco = 4
    dim_g = potencial_A.shape[1]
    F = np.zeros((dim_espaco, dim_espaco, dim_g), dtype=complex)

    for mu in range(dim_espaco):
        for nu in range(dim_espaco):
            # termo abeliano
            F[mu, nu, :] = (
                derivadas_parciais[mu, nu, :] - derivadas_parciais[nu, mu, :]
            )
            # termo não-abeliano −i g f^{abc} A_μ^b A_ν^c
            for a in range(dim_g):
                for b in range(dim_g):
                    for c in range(dim_g):
                        F[mu, nu, a] -= (
                            1j
                            * constante_acoplamento
                            * constantes_estrutura[b, c, a]
                            * potencial_A[mu, b]
                            * potencial_A[nu, c]
                        )
    return F


def exemplo_heisenberg_oscilador() -> Dict:
    """
    Álgebra de Heisenberg e oscilador harmônico
    (Cap. 1 §1.6.3 e Cap. 10 §10.2).
    """
    h3 = algebra_heisenberg()
    return {
        "algebra": str(h3),
        "relacoes": "[x,y]=z , [x,z]=0 , [y,z]=0",
        "aplicacao": "Quantização canônica do oscilador harmônico / princípio de incerteza",
        "eh_nilpotente": True,  # série central descendente termina (Cap. 4)
        "referencia": "Wilcke 2026, Cap. 1 §1.6.3 e Cap. 10 §10.2",
    }


def composicao_spins_meio() -> Dict:
    """
    Composição de dois spins 1/2 via Clebsch-Gordan
    (Cap. 10 §10.4 e Exercício 10.9).
    """
    cg = coeficientes_clebsch_gordan_su2(0.5, 0.5)
    return {
        "decomposicao": "1/2 ⊗ 1/2 = 0 ⊕ 1",
        "singleto": "estado singleto (antissimétrico) – singlete de spin",
        "tripleto": "estados tripleto (simétricos) – S=1",
        "coeficientes_exemplo": {
            "⟨↑↑|1,1⟩": cg.get((0.5, 0.5, 1.0, 1.0), None),
            "⟨↑↓|1,0⟩": cg.get((0.5, -0.5, 1.0, 0.0), None),
            "⟨↑↓|0,0⟩": cg.get((0.5, -0.5, 0.0, 0.0), None),
        },
        "referencia": "Wilcke 2026, Cap. 10 §10.4 eq. (10.13)",
    }
