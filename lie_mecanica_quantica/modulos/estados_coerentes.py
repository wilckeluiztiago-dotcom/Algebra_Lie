"""
Modelo 7 — Estados coerentes de Glauber
Autor: Luiz Tiago Wilcke
Fonte: Cap. 10 §10.2 (álgebra de Heisenberg / oscilador)
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple
from .oscilador_heisenberg import OsciladorHarmonico


def estado_coerente(
    alpha: complex,
    n_corte: int = 40,
) -> np.ndarray:
    """
    |α⟩ = e^{-|α|²/2} Σ (α^n / √n!) |n⟩
    Autovetor de a: a|α⟩ = α|α⟩
    """
    psi = np.zeros(n_corte + 1, dtype=complex)
    fator = np.exp(-0.5 * abs(alpha)**2)
    coef = fator
    psi[0] = coef
    for n in range(1, n_corte + 1):
        coef *= alpha / np.sqrt(n)
        psi[n] = coef
    # renormalizar por causa do corte
    norma = np.linalg.norm(psi)
    if norma > 0:
        psi /= norma
    return psi


def valor_esperado_coerente(
    alpha: complex,
    omega: float = 1.0,
    massa: float = 1.0,
    hbar: float = 1.0,
    n_corte: int = 40,
) -> Dict:
    """
    ⟨x⟩, ⟨p⟩, ⟨H⟩, ⟨N⟩ para um estado coerente.
    Valores teóricos:
      ⟨x⟩ = √(2ħ/(mω)) Re(α)
      ⟨p⟩ = √(2ħ m ω) Im(α)
      ⟨N⟩ = |α|²
      ⟨H⟩ = ħω (|α|² + 1/2)
    """
    osc = OsciladorHarmonico(omega, massa, hbar, n_corte)
    psi = estado_coerente(alpha, n_corte)
    # truncar operadores se necessário
    dim = min(len(psi), osc.dim)
    psi = psi[:dim]
    x = osc.x[:dim, :dim]
    p = osc.p[:dim, :dim]
    N = osc.N[:dim, :dim]
    H = osc.H[:dim, :dim]

    def ve(op):
        return float(np.real(psi.conj() @ op @ psi))

    return {
        "alpha": alpha,
        "⟨x⟩": ve(x),
        "⟨p⟩": ve(p),
        "⟨N⟩": ve(N),
        "⟨H⟩": ve(H),
        "⟨x⟩_teorico": np.sqrt(2 * hbar / (massa * omega)) * np.real(alpha),
        "⟨p⟩_teorico": np.sqrt(2 * hbar * massa * omega) * np.imag(alpha),
        "⟨N⟩_teorico": abs(alpha)**2,
        "⟨H⟩_teorico": hbar * omega * (abs(alpha)**2 + 0.5),
        "referencia": "Wilcke 2026, Cap. 10 §10.2",
    }


def sobreposicao_coerentes(alpha: complex, beta: complex, n_corte: int = 40) -> complex:
    """Produto escalar ⟨α|β⟩ = exp(−(|α|²+|β|²)/2 + α* β)."""
    return np.exp(-0.5 * (abs(alpha)**2 + abs(beta)**2) + np.conj(alpha) * beta)
