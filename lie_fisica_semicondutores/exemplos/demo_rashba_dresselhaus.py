#!/usr/bin/env python3
"""
Demonstração completa: Acoplamento spin-órbita Rashba + Dresselhaus.
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6 e Cap. 10 §10.3 do tratado
"""

import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.spin_orbita import (
    HamiltonianosSpinOrbita,
    espectro_rashba,
    espectro_rashba_dresselhaus,
    angulo_precessao_spin,
)
from modulos.grupos_classicos import algebra_su2


def main():
    print("=" * 70)
    print("  SPIN-ÓRBITA RASHBA + DRESSELHAUS — su(2) em 2DEG")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    # Parâmetros típicos de InAs ou InGaAs 2DEG
    m_ef = 0.023 * 9.109e-31
    alpha_R = 1.0e-11      # eV·m  (valor típico ~ 10^{-11} eV m)
    beta_D = 0.3e-11
    # converter para Joules·m
    alpha_J = alpha_R * 1.602e-19
    beta_J = beta_D * 1.602e-19

    print("\n[1] Parâmetros do 2DEG")
    print(f"    m* = 0.023 m_e")
    print(f"    α_R = {alpha_R:.2e} eV·m")
    print(f"    β_D = {beta_D:.2e} eV·m")

    modelo = HamiltonianosSpinOrbita(m_ef, alpha_J, beta_J)

    print("\n[2] Espectro em pontos selecionados do espaço k")
    print(f"    {'kx (1/nm)':>12} {'ky (1/nm)':>12} {'E- (meV)':>12} {'E+ (meV)':>12}")
    print("    " + "-" * 52)
    for kx_nm, ky_nm in [(0, 0), (0.1, 0), (0, 0.1), (0.1, 0.1), (0.2, 0)]:
        kx = kx_nm * 1e9
        ky = ky_nm * 1e9
        Em, Ep = modelo.espectro(kx, ky)
        print(f"    {kx_nm:12.2f} {ky_nm:12.2f} {Em/1.602e-22:12.4f} {Ep/1.602e-22:12.4f}")

    print("\n[3] Vetor de precessão Ω(k) e polarização de spin")
    kx, ky = 0.15e9, 0.05e9
    Omega = modelo.vetor_precessao(kx, ky)
    pol = modelo.polarizacao_spin(kx, ky)
    print(f"    Ω(k) = ({Omega[0]:.3e}, {Omega[1]:.3e}, {Omega[2]:.3e}) J")
    print(f"    polarização spin (estado inferior) = {pol}")

    print("\n[4] Ângulo de precessão (transistor de spin Datta-Das)")
    L = 1e-6  # 1 µm
    theta = angulo_precessao_spin(0.1e9, 0.0, alpha_J, beta_J, L)
    print(f"    Comprimento do canal L = 1 µm")
    print(f"    Ângulo de precessão θ ≈ {theta:.4f} rad  ({np.degrees(theta):.2f}°)")

    print("\n[5] Estrutura algébrica su(2)")
    su2 = algebra_su2()
    print(f"    {su2}")
    print(f"    Semissimples? {su2.eh_semissimples()}")
    print(f"    Constantes de estrutura ε_ijk (Jacobi satisfeita por construção)")

    print("\n[6] Espectro analítico puro de Rashba (β=0)")
    for k_nm in [0.0, 0.05, 0.10, 0.20]:
        k = k_nm * 1e9
        Em, Ep = espectro_rashba(k, m_ef, alpha_J)
        print(f"    k = {k_nm:.2f} nm⁻¹ → E± = {Em/1.602e-22:.4f} / {Ep/1.602e-22:.4f} meV")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 1 §1.6 e Cap. 10 §10.3")
    print("=" * 70)


if __name__ == "__main__":
    main()
