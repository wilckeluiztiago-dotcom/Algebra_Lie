"""
Módulo: representacoes
Representações irredutíveis, Clebsch-Gordan, pesos e fórmula de Weyl.

Referências no tratado:
- Capítulo 3: Representações de Álgebras de Lie
- §3.7 Decomposição Tensorial e o Teorema de Clebsch-Gordan
- Capítulo 9: Espaços de Pesos e Representações Irredutíveis
- §9.5 Fórmula de Dimensão de Weyl
- Cap. 10 §10.4 Produto de Clebsch-Gordan e Composição de Spins
"""

from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple, Optional
from math import factorial


def coeficientes_clebsch_gordan_su2(
    j1: float, j2: float
) -> Dict[Tuple[float, float, float, float], float]:
    """
    Coeficientes de Clebsch-Gordan ⟨j1 m1; j2 m2 | J M⟩ para su(2).

    Decomposição (tratado eq. 10.13):
    V^{j1} ⊗ V^{j2} = ⊕_{J=|j1-j2|}^{j1+j2} V^J

    Implementação analítica para os casos mais usados em física
    (spin 1/2 ⊗ 1/2, 1/2 ⊗ 1, etc.).
    """
    cg: Dict[Tuple[float, float, float, float], float] = {}

    # Caso canônico: 1/2 ⊗ 1/2 = 0 ⊕ 1  (singleto + tripleto)
    if abs(j1 - 0.5) < 1e-9 and abs(j2 - 0.5) < 1e-9:
        # singleto J=0
        cg[(0.5, 0.5, 0.0, 0.0)] = 1 / np.sqrt(2)
        cg[(0.5, -0.5, 0.0, 0.0)] = -1 / np.sqrt(2)
        # (m1,m2) permutados com sinal
        cg[(-0.5, 0.5, 0.0, 0.0)] = -1 / np.sqrt(2)  # já coberto pelo antissimétrico
        # tripleto J=1
        cg[(0.5, 0.5, 1.0, 1.0)] = 1.0
        cg[(-0.5, -0.5, 1.0, -1.0)] = 1.0
        cg[(0.5, -0.5, 1.0, 0.0)] = 1 / np.sqrt(2)
        cg[(-0.5, 0.5, 1.0, 0.0)] = 1 / np.sqrt(2)
        return cg

    # Caso genérico via fórmula de Racah (simplificada para j2=1/2)
    if abs(j2 - 0.5) < 1e-9:
        for m1 in np.arange(-j1, j1 + 0.1, 1.0):
            # J = j1 + 1/2
            J = j1 + 0.5
            M = m1 + 0.5
            if abs(M) <= J + 1e-9:
                cg[(m1, 0.5, J, M)] = np.sqrt((j1 + m1 + 1) / (2 * j1 + 1))
            M = m1 - 0.5
            if abs(M) <= J + 1e-9:
                cg[(m1, -0.5, J, M)] = np.sqrt((j1 - m1 + 1) / (2 * j1 + 1))
            # J = j1 - 1/2
            if j1 > 0:
                J = j1 - 0.5
                M = m1 + 0.5
                if abs(M) <= J + 1e-9:
                    cg[(m1, 0.5, J, M)] = -np.sqrt((j1 - m1) / (2 * j1 + 1))
                M = m1 - 0.5
                if abs(M) <= J + 1e-9:
                    cg[(m1, -0.5, J, M)] = np.sqrt((j1 + m1) / (2 * j1 + 1))
        return cg

    raise NotImplementedError(
        "Clebsch-Gordan genérico completo ainda não implementado; "
        "use os casos 1/2⊗1/2 ou j⊗1/2 (suficientes para o Cap. 10)."
    )


def pesos_representacao_su3(
    p: int, q: int
) -> List[Tuple[int, int]]:
    """
    Pesos (n1, n2) da representação irredutível (p,q) de su(3) = A₂
    (notação de Dynkin – Cap. 9 §9.2 e Fig. 9.1–9.2).

    Exemplo: (1,0) = fundamental 3
             (0,1) = antifundamental 3̄
             (1,1) = adjunta 8
             (3,0) = decupleto 10
    """
    pesos = []
    for i in range(p + q + 1):
        for j in range(p + q - i + 1):
            n1 = p - i + j
            n2 = q + i - 2 * j
            # filtro da regra de multiplicidade zero fora do hexágono/triângulo
            if n1 >= 0 and n2 >= 0 and n1 + n2 <= p + q:
                pesos.append((n1, n2))
    # versão simplificada: retorna o contorno principal
    # (implementação completa de multiplicidades exige fórmula de Freudenthal)
    return list(set(pesos))


def formula_dimensao_weyl(tipo: str, pesos_maximos: List[int]) -> int:
    """
    Fórmula de dimensão de Weyl (Cap. 9 §9.5):

    dim L(λ) = ∏_{α>0}  (⟨λ+ρ, α⟩) / (⟨ρ, α⟩)

    onde ρ = meia-soma das raízes positivas.

    Implementações explícitas para A_n e casos usados em física.
    """
    tipo = tipo.upper()
    if tipo == "A1":  # su(2)
        # λ = 2j → dim = 2j+1
        j = pesos_maximos[0] / 2
        return int(2 * j + 1)

    if tipo == "A2":  # su(3)
        # λ = p ω1 + q ω2
        p, q = pesos_maximos
        return int((p + 1) * (q + 1) * (p + q + 2) / 2)

    if tipo == "A3":  # su(4)
        a, b, c = pesos_maximos
        return int(
            (a + 1)
            * (b + 1)
            * (c + 1)
            * (a + b + 2)
            * (b + c + 2)
            * (a + b + c + 3)
            / 12
        )

    raise ValueError(f"Fórmula de Weyl para {tipo} ainda não codificada.")


def dimensao_irrep_su3(p: int, q: int) -> int:
    """Atalho conveniente – dim(p,q) de su(3)."""
    return formula_dimensao_weyl("A2", [p, q])
