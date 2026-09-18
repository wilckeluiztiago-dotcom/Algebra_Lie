"""
Modelo 9 — Teorema de Stone–von Neumann
Autor: Luiz Tiago Wilcke
Fonte: Cap. 10 §10.1 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple
from .oscilador_heisenberg import OsciladorHarmonico


def representacao_schrodinger(
    n_corte: int = 20,
    hbar: float = 1.0,
    massa: float = 1.0,
    omega: float = 1.0,
) -> Dict:
    """
    Representação de Schrödinger da álgebra de Heisenberg:
    x e p como operadores no espaço L², realizados na base de número
    do oscilador (equivalente unitária a qualquer outra irrep regular).

    O teorema de Stone–von Neumann (Cap. 10 §10.1) garante que
    qualquer representação irredutível “bonita” é unitariamente
    equivalente a esta.
    """
    osc = OsciladorHarmonico(omega, massa, hbar, n_corte)
    return {
        "x": osc.x,
        "p": osc.p,
        "a": osc.a,
        "a_dag": osc.a_dag,
        "dim": osc.dim,
        "enunciado": (
            "Toda representação irredutível da álgebra de Weyl-Heisenberg "
            "em espaço de Hilbert (com regularidade de Stone) é unitariamente "
            "equivalente à representação de Schrödinger."
        ),
        "referencia": "Wilcke 2026, Cap. 10 §10.1",
    }


def verificar_relacao_ccr(
    x: np.ndarray,
    p: np.ndarray,
    hbar: float = 1.0,
) -> Dict:
    """Verifica [x,p] ≈ i ħ na matriz truncada."""
    comut = x @ p - p @ x
    dim = x.shape[0]
    # bloco interior (a truncagem quebra a borda)
    b = slice(0, dim - 2)
    erro = np.max(np.abs(comut[b, b] - 1j * hbar * np.eye(dim - 2)))
    return {
        "max_erro_ccr": float(erro),
        "relacao": "[x, p] = i ħ",
        "referencia": "Wilcke 2026, Cap. 10 §10.1 e §10.2",
    }


def equivalencia_unitaria_exemplo() -> Dict:
    """
    Ilustração: a representação de momento (Fourier) é unitariamente
    equivalente à de posição via transformada de Fourier discreta.
    """
    N = 16
    # posição: diagonal
    x_pos = np.diag(np.linspace(-2, 2, N))
    # momento via DFT
    F = np.fft.fft(np.eye(N)) / np.sqrt(N)
    p_mom = F @ np.diag(np.linspace(-np.pi, np.pi, N)) @ F.conj().T
    return {
        "nota": "Representação de posição e de momento são unitariamente equivalentes (Stone–von Neumann)",
        "dim": N,
        "referencia": "Wilcke 2026, Cap. 10 §10.1",
    }
