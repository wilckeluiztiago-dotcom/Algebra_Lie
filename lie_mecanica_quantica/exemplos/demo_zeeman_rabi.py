#!/usr/bin/env python3
"""
Demo: Zeeman + Rabi (sistemas de dois níveis)
Autor: Luiz Tiago Wilcke — Cap. 10 §10.3
"""
import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.efeito_zeeman import HamiltonianoZeeman, espectro_zeeman, hamiltoniano_zeeman_spin12
from modulos.sistema_dois_niveis import HamiltonianoRabi, oscilacoes_rabi, frequencia_rabi


def main():
    print("=" * 70)
    print("  MODELOS 5 e 6 — Zeeman e Rabi")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    print("\n[1] Efeito Zeeman (spin 1/2, B ao longo de z)")
    Bz = 1.0
    E = espectro_zeeman(Bz, g=2.0, mu_B=1.0, j=0.5)
    print(f"    B = {Bz}  →  E = {E}  (esperado ±1)")

    print("\n[2] Zeeman com campo arbitrário")
    H = hamiltoniano_zeeman_spin12(0.3, 0.4, 0.5, omega_L=1.0)
    evals = np.linalg.eigvalsh(H)
    print(f"    B = (0.3,0.4,0.5), |B|={np.linalg.norm([0.3,0.4,0.5]):.4f}")
    print(f"    Autovalores = {evals}")

    print("\n[3] Frequência de Rabi")
    Omega, Delta = 1.0, 0.5
    Omega_R = frequencia_rabi(Omega, Delta)
    print(f"    Ω={Omega}, Δ={Delta} → Ω_R = {Omega_R:.4f}")

    print("\n[4] Oscilações de Rabi (estado inicial |g⟩)")
    rabi = HamiltonianoRabi(omega_0=1.0, Omega=1.0, omega=1.0)
    for t in [0.0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi]:
        P = rabi.probabilidade_excitacao(t)
        print(f"    t = {t:.3f}  →  P_e = {P:.4f}")

    print("\n[5] Curva analítica P_e(t)")
    curva = oscilacoes_rabi(Omega=1.0, Delta=0.0, t_max=4 * np.pi, n_pontos=5)
    print(f"    Ω_R = {curva['Omega_R']:.4f}")
    print(f"    P nos pontos amostrados: {curva['P_excitacao']}")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 10 §10.3")
    print("=" * 70)


if __name__ == "__main__":
    main()
