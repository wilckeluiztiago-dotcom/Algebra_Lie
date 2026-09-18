#!/usr/bin/env python3
"""
Testes dos 10 modelos de Mecânica Quântica com Álgebras de Lie.
Autor: Luiz Tiago Wilcke
"""
import sys
sys.path.insert(0, "..")
import numpy as np

from modulos.oscilador_heisenberg import OsciladorHarmonico, espectro_oscilador
from modulos.momento_angular_su2 import MomentoAngular
from modulos.adicao_spins_cg import adicionar_spins, base_acoplada
from modulos.casimir_invariantes import casimir_su2, verificar_comutacao_casimir
from modulos.efeito_zeeman import espectro_zeeman
from modulos.sistema_dois_niveis import frequencia_rabi, HamiltonianoRabi
from modulos.estados_coerentes import estado_coerente, valor_esperado_coerente
from modulos.evolucao_magnus_mq import evolucao_magnus, fidelidade_unitaria, exemplo_campo_magnetico_pulsado
from modulos.stone_von_neumann import representacao_schrodinger, verificar_relacao_ccr
from modulos.hidrogenio_angular import HarmonicosEsfericos, tabela_degenerescencia
from modulos.grupos_classicos import algebra_su2, algebra_heisenberg
from modulos.algebra_lie_basica import verificar_identidade_jacobi


def test_1_oscilador():
    E = espectro_oscilador(1.0, 3)
    assert np.allclose(E, [0.5, 1.5, 2.5, 3.5])
    osc = OsciladorHarmonico(n_corte=10)
    ver = osc.verificar_ccr()
    assert ver["max_|[a,a†] - 1|"] < 1e-10
    print("✓ 1 Oscilador / Heisenberg")


def test_2_momento_angular():
    ma = MomentoAngular(1.0)
    assert ma.dim == 3
    ver = ma.verificar_algebra()
    assert ver["max_|[Jx,Jy] - iħ Jz|"] < 1e-10
    print("✓ 2 Momento angular su(2)")


def test_3_cg():
    d = adicionar_spins(0.5, 0.5)
    assert d["dim_produto"] == d["dim_soma_irreps"] == 4
    base = base_acoplada(0.5, 0.5)
    # singleto normalizado
    s = base["estados_acoplados"][(0.0, 0.0)]
    assert abs(np.linalg.norm(s) - 1) < 1e-10
    print("✓ 3 Clebsch-Gordan")


def test_4_casimir():
    c = casimir_su2(1.0)
    assert abs(c["autovalor_C2"] - 2.0) < 1e-10
    com = verificar_comutacao_casimir(0.5)
    assert com["max_|[J²,Jx]|"] < 1e-10
    print("✓ 4 Casimir")


def test_5_zeeman():
    E = espectro_zeeman(1.0, g=2.0, mu_B=1.0, j=0.5)
    assert np.allclose(sorted(E), [-1.0, 1.0])
    print("✓ 5 Zeeman")


def test_6_rabi():
    assert abs(frequencia_rabi(3.0, 4.0) - 5.0) < 1e-10
    r = HamiltonianoRabi(1.0, 1.0, 1.0)
    P0 = r.probabilidade_excitacao(0.0)
    assert P0 < 1e-10
    print("✓ 6 Rabi")


def test_7_coerentes():
    psi = estado_coerente(1.0, n_corte=20)
    assert abs(np.linalg.norm(psi) - 1) < 1e-8
    ve = valor_esperado_coerente(1.0, n_corte=25)
    assert abs(ve["⟨N⟩"] - 1.0) < 0.05
    print("✓ 7 Estados coerentes")


def test_8_magnus():
    H = exemplo_campo_magnetico_pulsado(0.5, 0.3)
    U = evolucao_magnus(H, -1.0, 1.0, n_passos=40)
    assert fidelidade_unitaria(U) < 1e-10
    print("✓ 8 Magnus")


def test_9_stone():
    rep = representacao_schrodinger(n_corte=12)
    ccr = verificar_relacao_ccr(rep["x"], rep["p"])
    assert ccr["max_erro_ccr"] < 1e-8
    print("✓ 9 Stone–von Neumann")


def test_10_hidrogenio():
    he = HarmonicosEsfericos()
    esp = he.espectro_L2_Lz(2)
    assert abs(esp["autovalor_L2"] - 6.0) < 1e-10
    tab = tabela_degenerescencia(2)
    assert tab["degenerescencia_total_n"][1] == 1
    assert tab["degenerescencia_total_n"][2] == 4
    print("✓ 10 Hidrogênio angular")


def test_algebras():
    assert verificar_identidade_jacobi(algebra_su2().constantes_estrutura)
    assert verificar_identidade_jacobi(algebra_heisenberg().constantes_estrutura)
    print("✓ Álgebras su(2) e Heisenberg")


if __name__ == "__main__":
    test_1_oscilador()
    test_2_momento_angular()
    test_3_cg()
    test_4_casimir()
    test_5_zeeman()
    test_6_rabi()
    test_7_coerentes()
    test_8_magnus()
    test_9_stone()
    test_10_hidrogenio()
    test_algebras()
    print("\nTodos os 10 modelos passaram nos testes.")
