#!/usr/bin/env python3
"""
Demonstração: su(2) / momento angular / isospin
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6 e Cap. 10 §10.3 do tratado
"""

import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.grupos_classicos import algebra_su2, geradores_pauli
from modulos.representacoes import coeficientes_clebsch_gordan_su2


def main():
    print("=" * 60)
    print("  su(2) — Momento Angular e Isospin")
    print("  (Wilcke 2026, Cap. 1 §1.6.1–1.6.2 e Cap. 10 §10.3)")
    print("=" * 60)

    su2 = algebra_su2()
    print(f"\nÁlgebra: {su2}")
    print(f"Base: {su2.base_nomes}")
    print(f"Semissimples? {su2.eh_semissimples()}")

    # Constantes de estrutura ε_ijk
    print("\nConstantes de estrutura não-nulas (ε_ijk):")
    c = su2.constantes_estrutura
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if abs(c[i, j, k]) > 1e-10:
                    print(f"  c^{k}_{i}{j} = {c[i,j,k]:.0f}")

    # Matrizes de Pauli
    print("\nGeradores de Pauli (base hermitiana):")
    for i, sigma in enumerate(geradores_pauli()):
        print(f"  σ_{i+1} =\n{sigma}")

    # Clebsch-Gordan 1/2 ⊗ 1/2
    print("\nClebsch-Gordan 1/2 ⊗ 1/2 = 0 ⊕ 1:")
    cg = coeficientes_clebsch_gordan_su2(0.5, 0.5)
    for chave, val in sorted(cg.items()):
        m1, m2, J, M = chave
        print(f"  ⟨{m1:+.1f},{m2:+.1f}|{J:.0f},{M:+.1f}⟩ = {val:.4f}")

    print("\nReferência: tratado de Luiz Tiago Wilcke, 2026.")


if __name__ == "__main__":
    main()
