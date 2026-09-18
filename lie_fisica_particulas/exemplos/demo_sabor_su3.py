#!/usr/bin/env python3
"""
Demonstração: SU(3) de sabor – octeto e decupleto
Autor: Luiz Tiago Wilcke
Fonte: Cap. 9 e Cap. 10 §10.6.1
"""

import sys
sys.path.insert(0, "..")
from modulos.grupos_classicos import algebra_su3, geradores_gell_mann
from modulos.representacoes import dimensao_irrep_su3, formula_dimensao_weyl
from modulos.sistemas_raizes import sistema_raizes_a2, matriz_de_cartan, diagrama_dynkin_simples
from modulos.fisica_particulas import decomposicao_mesons, decomposicao_barions


def main():
    print("=" * 65)
    print("  SU(3) de sabor – Mésons e Bárions")
    print("  (Wilcke 2026, Cap. 9 e Cap. 10 §10.6.1)")
    print("=" * 65)

    su3 = algebra_su3()
    print(f"\nÁlgebra: {su3}")
    print(f"dim su(3) = {su3.dimensao} → 8 glúons (ou 8 mésons do octeto)")

    print("\nMatriz de Cartan de A₂:")
    print(matriz_de_cartan("A2"))
    print("Diagrama de Dynkin:", diagrama_dynkin_simples("A2"))

    print("\nDimensões de irreps (fórmula de Weyl):")
    for p, q, nome in [(1, 0, "3"), (0, 1, "3̄"), (1, 1, "8"), (3, 0, "10"), (2, 2, "27")]:
        d = dimensao_irrep_su3(p, q)
        print(f"  ({p},{q}) = {nome:>3}  →  dim = {d}")

    print("\n--- Mésons ---")
    mes = decomposicao_mesons()
    print(f"{mes['produto']} = {mes['decomposicao']}")
    print("Partículas:", ", ".join(mes["particulas_octeto"]))

    print("\n--- Bárions ---")
    bar = decomposicao_barions()
    print(f"{bar['produto']} = {bar['decomposicao']}")
    print(bar["destaque"])

    print("\nSistema de raízes A₂ (6 raízes):")
    phi = sistema_raizes_a2()
    print(f"  posto = {phi.posto}, |Φ| = {phi.numero_raizes}")

    print("\nReferência completa: tratado de Luiz Tiago Wilcke, 2026.")


if __name__ == "__main__":
    main()
