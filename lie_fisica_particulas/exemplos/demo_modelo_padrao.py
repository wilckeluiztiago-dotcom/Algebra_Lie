#!/usr/bin/env python3
"""
Demonstração: Modelo Padrão a partir das álgebras de Lie do tratado.
Autor: Luiz Tiago Wilcke
Fonte: Cap. 10 do tratado Álgebras de Lie (2026)
"""

import sys
sys.path.insert(0, "..")

from modulos.fisica_particulas import (
    modelo_padrao,
    decomposicao_mesons,
    decomposicao_barions,
    grupo_gut_su5,
    cadeia_excepcional_e8,
    exemplo_heisenberg_oscilador,
    composicao_spins_meio,
)
from modulos.grupos_classicos import algebra_su3, algebra_su2
from modulos.algebra_lie_basica import verificar_identidade_jacobi


def main():
    print("=" * 70)
    print("  ÁLGEBRAS DE LIE → FÍSICA DE PARTÍCULAS")
    print("  Autor do tratado e do código: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    # 1. Modelo Padrão
    print("\n[1] Modelo Padrão SU(3)×SU(2)×U(1)  (Cap. 10 §10.6)")
    sm = modelo_padrao()
    for k, v in sm.items():
        print(f"    {k}: {v}")

    # 2. Verificação da identidade de Jacobi para su(3)
    print("\n[2] Verificação da Identidade de Jacobi em su(3)")
    su3 = algebra_su3()
    ok = verificar_identidade_jacobi(su3.constantes_estrutura)
    print(f"    Identidade de Jacobi satisfeita? {ok}")
    print(f"    Álgebra semissimples (det κ ≠ 0)? {su3.eh_semissimples()}")

    # 3. Mésons e Bárions
    print("\n[3] Decomposição de Mésons 3⊗3̄ = 1⊕8  (eq. 10.28)")
    mes = decomposicao_mesons()
    print(f"    {mes['produto']} = {mes['decomposicao']}")
    print(f"    Dimensões: {mes['dimensoes']}")
    print(f"    Partículas do octeto: {', '.join(mes['particulas_octeto'])}")

    print("\n[4] Decomposição de Bárions 3⊗3⊗3 = 1⊕8⊕8⊕10  (eq. 10.29)")
    bar = decomposicao_barions()
    print(f"    {bar['produto']} = {bar['decomposicao']}")
    print(f"    Destaque: {bar['destaque']}")

    # 5. GUT SU(5)
    print("\n[5] GUT SU(5) de Georgi-Glashow  (Cap. 10 §10.7)")
    gut = grupo_gut_su5()
    print(f"    dim SU(5) = {gut['dimensao']}")
    print(f"    Conteúdo: {gut['conteudo']}")
    print(f"    Decaimento do próton: {gut['decaimento_proton']}")

    # 6. Cadeia excepcional
    print("\n[6] Cadeia excepcional E₈ ⊃ … ⊃ SM  (Cap. 8 e 10)")
    cad = cadeia_excepcional_e8()
    print(f"    {cad['cadeia']}")
    print(f"    Dimensões: {cad['dimensoes']}")

    # 7. Heisenberg
    print("\n[7] Álgebra de Heisenberg (Cap. 1 §1.6.3 + Cap. 10 §10.2)")
    h = exemplo_heisenberg_oscilador()
    print(f"    {h['algebra']}")
    print(f"    Relações: {h['relacoes']}")

    # 8. Spins
    print("\n[8] Composição de spins 1/2 ⊗ 1/2  (Cap. 10 §10.4)")
    spins = composicao_spins_meio()
    print(f"    {spins['decomposicao']}")
    print(f"    Coeficientes CG exemplo: {spins['coeficientes_exemplo']}")

    print("\n" + "=" * 70)
    print("  Todas as referências estão no tratado de Luiz Tiago Wilcke (2026).")
    print("=" * 70)


if __name__ == "__main__":
    main()
