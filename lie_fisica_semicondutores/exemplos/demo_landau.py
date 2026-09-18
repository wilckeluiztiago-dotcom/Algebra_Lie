#!/usr/bin/env python3
"""
Demonstração completa: Níveis de Landau a partir da álgebra de Heisenberg.
Autor: Luiz Tiago Wilcke
Fonte: Cap. 1 §1.6.3 e Cap. 10 §10.2 do tratado Álgebras de Lie (2026)
"""

import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.niveis_landau import (
    espectro_landau,
    HamiltonianosLandau,
    operadores_ladder_landau,
)


def main():
    print("=" * 70)
    print("  NÍVEIS DE LANDAU — Álgebra de Heisenberg em Semicondutores")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    # Parâmetros típicos de GaAs 2DEG
    B = 5.0          # Tesla
    m_ef = 0.067 * 9.109e-31  # massa efetiva GaAs
    n_max = 5

    print("\n[1] Parâmetros do 2DEG (GaAs)")
    print(f"    B = {B} T")
    print(f"    m* = {m_ef:.3e} kg  (0.067 m_e)")

    omega_c, ell_B = operadores_ladder_landau(B, m_ef)
    print(f"\n[2] Frequência ciclotrônica e comprimento magnético")
    print(f"    ω_c = {omega_c:.4e} rad/s")
    print(f"    ℓ_B = {ell_B*1e9:.3f} nm")
    print(f"    ħω_c = {1.054e-34 * omega_c / 1.602e-19 * 1e3:.3f} meV")

    print("\n[3] Espectro de Landau (n = 0…5) com efeito Zeeman")
    espec = espectro_landau(n_max, B, m_ef, g_fator=0.44)  # g* GaAs ≈ 0.44
    print(f"    {'n':>3} {'spin':>5} {'E (meV)':>12}")
    print("    " + "-" * 25)
    for niv in espec["niveis"]:
        print(f"    {niv['n']:3d} {niv['spin']:+5d} {niv['E_total_eV']*1e3:12.4f}")

    print("\n[4] Verificação da origem Heisenberg")
    H_land = HamiltonianosLandau(n_corte=4, B=B, m_efetiva=m_ef)
    verif = H_land.verificar_algebra_heisenberg()
    print(f"    Álgebra: {verif['algebra']}")
    print(f"    [x,y] → {verif['comutador_[x,y]']}")
    print(f"    Relação: {verif['relacao_landau']}")

    print("\n[5] Densidade de estados por nível")
    area = 1e-12  # 1 µm²
    deg = H_land.densidade_estados_landau(area)
    print(f"    Área = 1 µm² → degenerescência ≈ {deg:.1f} estados/nível")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 1 §1.6.3 e Cap. 10 §10.2")
    print("=" * 70)


if __name__ == "__main__":
    main()
