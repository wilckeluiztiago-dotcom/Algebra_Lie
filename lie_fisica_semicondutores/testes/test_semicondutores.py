#!/usr/bin/env python3
"""
Testes unitários completos — Física de Semicondutores com Álgebras de Lie.
Autor: Luiz Tiago Wilcke
"""

import sys
sys.path.insert(0, "..")
import numpy as np

from modulos.niveis_landau import espectro_landau, HamiltonianosLandau, operadores_ladder_landau
from modulos.spin_orbita import HamiltonianosSpinOrbita, espectro_rashba
from modulos.evolucao_magnus import (
    exemplo_pulso_eletrico,
    evolucao_unitaria_magnus,
    fidelidade_unitaria,
    comparacao_euler_vs_magnus,
)
from modulos.excitons_cg import estados_exciton_cg, hamiltoniano_estrutura_fina
from modulos.vales_su2 import HamiltonianosVale, espectro_grafeno_massivo
from modulos.grupos_classicos import algebra_su2, algebra_heisenberg
from modulos.algebra_lie_basica import verificar_identidade_jacobi


def test_heisenberg_landau():
    h3 = algebra_heisenberg()
    assert verificar_identidade_jacobi(h3.constantes_estrutura)
    omega, ell = operadores_ladder_landau(1.0, 9.1e-31)
    assert omega > 0 and ell > 0
    espec = espectro_landau(2, 1.0, 9.1e-31)
    assert len(espec["niveis"]) == 6  # 3 n × 2 spins
    print("✓ Heisenberg + níveis de Landau OK")


def test_rashba():
    m = 0.067 * 9.1e-31
    alpha = 1e-11 * 1.602e-19
    modelo = HamiltonianosSpinOrbita(m, alpha, 0.0)
    Em, Ep = modelo.espectro(1e8, 0.0)
    Em2, Ep2 = espectro_rashba(1e8, m, alpha)
    assert abs(Em - Em2) < 1e-25
    assert abs(Ep - Ep2) < 1e-25
    # gap de Rashba = 2 α k
    assert abs((Ep - Em) - 2 * alpha * 1e8) < 1e-25
    print("✓ Rashba / su(2) OK")


def test_magnus_unitario():
    H = exemplo_pulso_eletrico(1e-22, 1e-12, 9.1e-31, alpha_R=1e-23)
    U = evolucao_unitaria_magnus(H, -3e-12, 3e-12, n_passos=32)
    fid = fidelidade_unitaria(U)
    assert fid < 1e-10
    psi0 = np.array([1.0, 0.0], dtype=complex)
    comp = comparacao_euler_vs_magnus(H, psi0, -1e-12, 1e-12, n_passos=20)
    assert comp["erro_norma_lie"] < 1e-10
    print("✓ Magnus unitário e Lie-Euler OK")


def test_exciton_cg():
    est = estados_exciton_cg(0.5, 1.5)
    assert "1 ⊕ 2" in est["decomposicao"] or "1/2 ⊗ 3/2" in est["decomposicao"]
    H = hamiltoniano_estrutura_fina(0.001, -0.002, 0.0005)
    assert H.shape == (8, 8)
    assert np.allclose(H, H.conj().T)  # hermitiano
    print("✓ Éxciton Clebsch-Gordan OK")


def test_vales():
    vale = HamiltonianosVale(1e6, 0.1 * 1.602e-19)
    Em, Ep = vale.espectro(0.0, 0.0)
    assert abs(Em + 0.1 * 1.602e-19) < 1e-25
    assert abs(Ep - 0.1 * 1.602e-19) < 1e-25
    Em2, Ep2 = espectro_grafeno_massivo(1e8, 1e6, 0.1 * 1.602e-19)
    assert Em2 < 0 < Ep2
    print("✓ Vales su(2) / grafeno OK")


def test_su2_semissimples():
    su2 = algebra_su2()
    assert su2.eh_semissimples()
    assert verificar_identidade_jacobi(su2.constantes_estrutura)
    print("✓ su(2) semissimples OK")


if __name__ == "__main__":
    test_heisenberg_landau()
    test_rashba()
    test_magnus_unitario()
    test_exciton_cg()
    test_vales()
    test_su2_semissimples()
    print("\nTodos os testes de semicondutores passaram.")
