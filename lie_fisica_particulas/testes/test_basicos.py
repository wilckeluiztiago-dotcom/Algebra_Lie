#!/usr/bin/env python3
"""
Testes unitários básicos dos módulos de Álgebras de Lie.
Autor: Luiz Tiago Wilcke
"""

import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.grupos_classicos import algebra_su2, algebra_su3, algebra_heisenberg
from modulos.algebra_lie_basica import verificar_identidade_jacobi
from modulos.representacoes import dimensao_irrep_su3, formula_dimensao_weyl
from modulos.fisica_particulas import modelo_padrao, decomposicao_mesons, decomposicao_barions


def test_jacobi_su2():
    su2 = algebra_su2()
    assert verificar_identidade_jacobi(su2.constantes_estrutura)
    print("✓ Jacobi su(2) OK")


def test_jacobi_su3():
    su3 = algebra_su3()
    assert verificar_identidade_jacobi(su3.constantes_estrutura)
    print("✓ Jacobi su(3) OK")


def test_heisenberg_nilpotente():
    h = algebra_heisenberg()
    # [x,[x,y]] = [x,z] = 0 → nilpotente
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = h.colchete_lie(x, y)
    assert np.allclose(z, [0, 0, 1])
    assert np.allclose(h.colchete_lie(x, z), 0)
    print("✓ Heisenberg [x,y]=z e central OK")


def test_dimensoes_su3():
    assert dimensao_irrep_su3(1, 0) == 3
    assert dimensao_irrep_su3(1, 1) == 8
    assert dimensao_irrep_su3(3, 0) == 10
    assert formula_dimensao_weyl("A1", [1]) == 2  # j=1/2
    print("✓ Dimensões de Weyl OK")


def test_modelo_padrao():
    sm = modelo_padrao()
    assert sm["total_bosons_gauge"] == 12
    assert sm["gluons"] == 8
    print("✓ Modelo Padrão 12 bósons OK")


def test_decomposicoes():
    mes = decomposicao_mesons()
    assert mes["dimensoes"]["8"] == 8
    bar = decomposicao_barions()
    assert bar["dimensoes"]["10"] == 10
    print("✓ Decomposições 3⊗3̄ e 3⊗3⊗3 OK")


if __name__ == "__main__":
    test_jacobi_su2()
    test_jacobi_su3()
    test_heisenberg_nilpotente()
    test_dimensoes_su3()
    test_modelo_padrao()
    test_decomposicoes()
    print("\nTodos os testes passaram.")
