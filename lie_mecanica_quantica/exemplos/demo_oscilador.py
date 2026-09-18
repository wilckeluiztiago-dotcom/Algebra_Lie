#!/usr/bin/env python3
"""
Demo: Oscilador harmônico + Heisenberg + estados coerentes
Autor: Luiz Tiago Wilcke — Cap. 1 §1.6.3 e Cap. 10 §10.2
"""
import sys
sys.path.insert(0, "..")
import numpy as np
from modulos.oscilador_heisenberg import OsciladorHarmonico, verificar_heisenberg_algebra, espectro_oscilador
from modulos.estados_coerentes import estado_coerente, valor_esperado_coerente
from modulos.stone_von_neumann import representacao_schrodinger, verificar_relacao_ccr


def main():
    print("=" * 70)
    print("  MODELOS 1, 7 e 9 — Oscilador, Coerentes, Stone–von Neumann")
    print("  Autor: Luiz Tiago Wilcke (2026)")
    print("=" * 70)

    print("\n[1] Álgebra de Heisenberg abstrata")
    h = verificar_heisenberg_algebra()
    print(f"    {h['algebra']}")
    print(f"    [x,y] = {h['[x,y]']}")
    print(f"    {h['relacao']}")

    print("\n[2] Espectro do oscilador (ħ=ω=1)")
    E = espectro_oscilador(1.0, 5)
    for n, e in enumerate(E):
        print(f"    E_{n} = {e:.4f}")

    osc = OsciladorHarmonico(omega=1.0, n_corte=20)
    ver = osc.verificar_ccr()
    print(f"\n[3] Verificação CCR (truncagem)")
    print(f"    max |[x,p]−iħ| ≈ {ver['max_|[x,p] - iħ|']:.2e}")
    print(f"    max |[a,a†]−1| ≈ {ver['max_|[a,a†] - 1|']:.2e}")

    print("\n[4] Estado coerente α = 1.5 + 0.5i")
    alpha = 1.5 + 0.5j
    ve = valor_esperado_coerente(alpha, n_corte=30)
    print(f"    ⟨N⟩ = {ve['⟨N⟩']:.4f}  (teórico {ve['⟨N⟩_teorico']:.4f})")
    print(f"    ⟨H⟩ = {ve['⟨H⟩']:.4f}  (teórico {ve['⟨H⟩_teorico']:.4f})")
    print(f"    ⟨x⟩ = {ve['⟨x⟩']:.4f}  (teórico {ve['⟨x⟩_teorico']:.4f})")

    print("\n[5] Stone–von Neumann")
    rep = representacao_schrodinger(n_corte=15)
    print(f"    {rep['enunciado'][:80]}...")
    ccr = verificar_relacao_ccr(rep["x"], rep["p"])
    print(f"    max erro CCR = {ccr['max_erro_ccr']:.2e}")

    print("\n" + "=" * 70)
    print("  Referência: Wilcke, Álgebras de Lie, Cap. 1 e Cap. 10")
    print("=" * 70)


if __name__ == "__main__":
    main()
