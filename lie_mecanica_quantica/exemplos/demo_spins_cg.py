#!/usr/bin/env python3
"""
Demo: Momento angular, Casimir, Clebsch-Gordan, Hidrogênio
Autor: Luiz Tiago Wilcke — Cap. 1, 3 e 10
"""
import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.momento_angular_su2 import MomentoAngular, autovalores_j2_jz
from modulos.adicao_spins_cg import adicionar_spins, base_acoplada
from modulos.casimir_invariantes import casimir_su2, verificar_comutacao_casimir
from modulos.hidrogenio_angular import HarmonicosEsfericos, tabela_degenerescencia


def main():
    print("=" * 70)
    print("  MODELOS 2, 3, 4 e 10 — su(2), CG, Casimir, Hidrogênio")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    print("\n[1] Momento angular j=1")
    ma = MomentoAngular(1.0)
    print(f"    dim = {ma.dim}")
    print(f"    autovalores Jz = {ma.autovalores_jz()}")
    print(f"    autovalor J² = {ma.autovalores_j2()[0]:.4f}")
    ver = ma.verificar_algebra()
    print(f"    max |[Jx,Jy]−iħJz| = {ver['max_|[Jx,Jy] - iħ Jz|']:.2e}")

    print("\n[2] Casimir")
    c = casimir_su2(1.0)
    print(f"    C₂ autovalor = {c['autovalor_C2']:.4f}")
    print(f"    desvio da identidade = {c['max_|C2 - j(j+1)ħ² 1|']:.2e}")
    com = verificar_comutacao_casimir(1.0)
    print(f"    [J²,Jx] max = {com['max_|[J²,Jx]|']:.2e}")

    print("\n[3] Adição 1/2 ⊗ 1/2")
    dec = adicionar_spins(0.5, 0.5)
    print(f"    {dec['decomposicao']}")
    base = base_acoplada(0.5, 0.5)
    print(f"    |1,0⟩ = {base['estados_acoplados'][(1.0, 0.0)]}")
    print(f"    |0,0⟩ = {base['estados_acoplados'][(0.0, 0.0)]}")

    print("\n[4] Parte angular do hidrogênio")
    he = HarmonicosEsfericos(l_max=2)
    for l in range(3):
        esp = he.espectro_L2_Lz(l)
        print(f"    l={l}: L²={esp['autovalor_L2']:.1f}, dim={esp['dim']}")
    tab = tabela_degenerescencia(3)
    print(f"    Degenerescência n=1,2,3: {tab['degenerescencia_total_n']}")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 1, 3 e 10")
    print("=" * 70)


if __name__ == "__main__":
    main()
