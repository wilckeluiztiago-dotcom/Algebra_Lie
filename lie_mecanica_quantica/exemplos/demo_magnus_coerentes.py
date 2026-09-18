#!/usr/bin/env python3
"""
Demo: Magnus + coerentes (evolução temporal)
Autor: Luiz Tiago Wilcke — Cap. 12 e Cap. 10
"""
import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.evolucao_magnus_mq import (
    exemplo_campo_magnetico_pulsado,
    evolucao_magnus,
    fidelidade_unitaria,
    comparacao_preservacao_norma,
)
from modulos.estados_coerentes import estado_coerente, sobreposicao_coerentes


def main():
    print("=" * 70)
    print("  MODELOS 7 e 8 — Coerentes + Magnus")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    print("\n[1] Sobreposição de estados coerentes")
    alpha, beta = 1.0 + 0.5j, 0.5 - 0.3j
    prod = sobreposicao_coerentes(alpha, beta)
    print(f"    ⟨α|β⟩ = {prod:.6f}")
    print(f"    |⟨α|β⟩|² = {abs(prod)**2:.6f}")

    print("\n[2] Evolução com Magnus (pulso Gaussiano de campo B)")
    H = exemplo_campo_magnetico_pulsado(B0=1.0, tau=0.5, direcao="x")
    U = evolucao_magnus(H, t0=-2.0, t1=2.0, n_passos=80)
    fid = fidelidade_unitaria(U)
    print(f"    Fidelidade unitária ‖U†U−1‖ = {fid:.2e}")

    psi0 = np.array([1.0, 0.0], dtype=complex)
    psi_f = U @ psi0
    print(f"    |ψ₀⟩ = |↑⟩  →  |ψ_f|² = {np.abs(psi_f)**2}")

    print("\n[3] Comparação Lie-Euler vs Euler ordinário (Fig. 12.1)")
    comp = comparacao_preservacao_norma(H, psi0, -1.0, 1.0, n_passos=40)
    print(f"    Erro norma Lie-Euler : {comp['erro_lie']:.2e}")
    print(f"    Erro norma Euler ord.: {comp['erro_euler']:.2e}")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 12 §12.2–12.4")
    print("=" * 70)


if __name__ == "__main__":
    main()
