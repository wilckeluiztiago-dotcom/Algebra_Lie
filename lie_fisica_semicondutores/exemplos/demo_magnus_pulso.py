#!/usr/bin/env python3
"""
Demonstração completa: Série de Magnus e integradores em grupos de Lie
para evolução temporal sob pulsos em semicondutores.
Autor: Luiz Tiago Wilcke
Fonte: Capítulo 12 do tratado Álgebras de Lie (2026)
"""

import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.evolucao_magnus import (
    exemplo_pulso_eletrico,
    evolucao_unitaria_magnus,
    integrador_lie_euler,
    comparacao_euler_vs_magnus,
    fidelidade_unitaria,
    serie_magnus_ordem2,
)


def main():
    print("=" * 70)
    print("  SÉRIE DE MAGNUS — Evolução unitária sob pulsos")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    # Pulso Gaussiano de Zeeman + termo de Rashba constante
    E0 = 0.5e-3 * 1.602e-19   # 0.5 meV em Joules
    tau = 1e-12               # 1 ps
    alpha_R = 0.1e-3 * 1.602e-19

    H = exemplo_pulso_eletrico(E0, tau, m_efetiva=0.067*9.1e-31, alpha_R=alpha_R)

    print("\n[1] Hamiltoniano de pulso")
    print("    H(t) = E0 exp(−t²/(2τ²)) σ_z  +  α_R σ_x")
    print(f"    E0 = 0.5 meV,  τ = 1 ps,  α_R = 0.1 meV")

    t0, t1 = -3 * tau, 3 * tau
    print(f"\n[2] Intervalo de integração: t ∈ [{t0*1e12:.1f}, {t1*1e12:.1f}] ps")

    # Magnus ordem 2
    print("\n[3] Gerador de Magnus (ordem 2)")
    Omega = serie_magnus_ordem2(H, t0, t1, n_passos=128)
    print(f"    ‖Ω‖_F = {np.linalg.norm(Omega, ord='fro'):.6e}")
    print(f"    Ω é anti-hermitiano? ‖Ω + Ω†‖ = {np.linalg.norm(Omega + Omega.conj().T):.2e}")

    U = evolucao_unitaria_magnus(H, t0, t1, n_passos=128)
    fid = fidelidade_unitaria(U)
    print(f"\n[4] Operador de evolução U = exp(Ω)")
    print(f"    Fidelidade unitária ‖U†U − 1‖_F = {fid:.2e}  (deve ser ~0)")

    # Estado inicial |↑⟩
    psi0 = np.array([1.0, 0.0], dtype=complex)
    psi_final = U @ psi0
    print(f"\n[5] Evolução do estado |↑⟩")
    print(f"    |ψ_final|² = {np.abs(psi_final)**2}")
    print(f"    Probabilidade de |↓⟩ = {np.abs(psi_final[1])**2:.6f}")

    # Comparação Euler vs Lie
    print("\n[6] Comparação integrador de Lie vs Euler ordinário (Fig. 12.1 do tratado)")
    comp = comparacao_euler_vs_magnus(H, psi0, t0, t1, n_passos=80)
    print(f"    Erro máximo de norma (Lie-Euler)  : {comp['erro_norma_lie']:.2e}")
    print(f"    Erro máximo de norma (Euler ord.) : {comp['erro_norma_euler']:.2e}")
    print("    → O integrador geométrico preserva a norma; o Euler ordinário drena.")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 12 §12.2 e §12.4")
    print("=" * 70)


if __name__ == "__main__":
    main()
