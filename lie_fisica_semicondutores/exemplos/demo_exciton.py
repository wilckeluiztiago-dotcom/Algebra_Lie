#!/usr/bin/env python3
"""
Demonstração completa: Estrutura fina de éxcitons via Clebsch-Gordan.
Autor: Luiz Tiago Wilcke
Fonte: Cap. 3 §3.7 e Cap. 10 §10.4 do tratado
"""

import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.excitons_cg import (
    estados_exciton_cg,
    hamiltoniano_estrutura_fina,
    espectro_estrutura_fina,
    intensidade_optica_relativa,
)
from modulos.vales_su2 import (
    HamiltonianosVale,
    espectro_grafeno_massivo,
    densidade_estados_2d_dirac,
)


def main():
    print("=" * 70)
    print("  ÉXCITONS E VALES — Clebsch-Gordan e su(2) de vale")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    # ---- Éxcitons ----
    print("\n[1] Acoplamento elétron (j=1/2) ⊗ buraco (j=3/2)")
    est = estados_exciton_cg(0.5, 1.5)
    print(f"    Decomposição: {est['decomposicao']}")
    print(f"    J possíveis: {est['J_possiveis']}")
    print(f"    Interpretação:")
    for k, v in est["interpretacao"].items():
        print(f"      {k}: {v}")

    print("\n[2] Alguns coeficientes de Clebsch-Gordan")
    for chave, val in list(est["coeficientes_CG"].items())[:6]:
        me, mh, J, M = chave
        print(f"    ⟨{me:+.1f},{mh:+.1f}|{J:.0f},{M:+.1f}⟩ = {val:.4f}")

    print("\n[3] Hamiltoniano de estrutura fina e espectro")
    delta0 = 0.5e-3   # eV
    delta1 = -1.0e-3
    delta2 = 0.2e-3
    espec = espectro_estrutura_fina(delta0, delta1, delta2)
    print(f"    δ₀ = {delta0*1e3:.2f} meV,  δ₁ = {delta1*1e3:.2f} meV,  δ₂ = {delta2*1e3:.2f} meV")
    print(f"    Energias (meV): {espec['energias']*1e3}")

    print("\n[4] Intensidades ópticas relativas (σ+)")
    intens = intensidade_optica_relativa(est, "sigma+")
    for estado, I in intens.items():
        print(f"    {estado}: I_rel = {I:.4f}")

    # ---- Vales (grafeno / TMD) ----
    print("\n" + "-" * 70)
    print("[5] Hamiltoniano de vale (grafeno massivo / TMD)")
    v_F = 1.0e6          # m/s
    Delta = 0.1 * 1.602e-19  # 0.1 eV em J
    vale = HamiltonianosVale(v_F, Delta)

    print(f"    v_F = {v_F:.1e} m/s,  Δ = 0.1 eV")
    for k_nm in [0.0, 0.05, 0.10, 0.20]:
        k = k_nm * 1e9
        Em, Ep = vale.espectro(k, 0.0, tau=+1)
        print(f"    k = {k_nm:.2f} nm⁻¹ → E± = {Em/1.602e-19:.4f} / {Ep/1.602e-19:.4f} eV")

    print("\n[6] Densidade de estados 2D Dirac")
    for E_eV in [0.05, 0.10, 0.20, 0.50]:
        D = densidade_estados_2d_dirac(E_eV * 1.602e-19, v_F, Delta)
        print(f"    D(E={E_eV:.2f} eV) = {D:.4e}  (estados / (J·m²))")

    print("\n[7] Matriz 4×4 (dois vales)")
    H4 = vale.hamiltoniano_quatro_componentes(0.1e9, 0.05e9)
    evals4 = np.linalg.eigvalsh(H4)
    print(f"    Autovalores (eV): {evals4 / 1.602e-19}")

    print("\n" + "=" * 70)
    print("  Referências: Wilcke 2026, Cap. 3 §3.7, Cap. 10 §10.4 e Cap. 1 §1.6")
    print("=" * 70)


if __name__ == "__main__":
    main()
