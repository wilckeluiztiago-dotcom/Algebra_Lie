"""
Modelo 8 — Evolução temporal via série de Magnus
Autor: Luiz Tiago Wilcke
Fonte: Cap. 12 §12.2–12.4 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Callable, Dict, Tuple
from scipy.linalg import expm


def serie_magnus_ordem2(
    H: Callable[[float], np.ndarray],
    t0: float,
    t1: float,
    n_passos: int = 64,
    hbar: float = 1.0,
) -> np.ndarray:
    """
    Ω ≈ Ω₁ + Ω₂  (Cap. 12 §12.2)

    Ω₁ = −(i/ħ) ∫ H(t) dt
    Ω₂ = (1/2) (−i/ħ)² ∫∫ [H(t1), H(t2)] dt1 dt2
    """
    dt = (t1 - t0) / n_passos
    tempos = np.linspace(t0 + dt / 2, t1 - dt / 2, n_passos)
    Omega1 = np.zeros_like(H(t0), dtype=complex)
    for t in tempos:
        Omega1 += -1j / hbar * H(t) * dt

    Omega2 = np.zeros_like(Omega1)
    for i, t1_ in enumerate(tempos):
        H1 = H(t1_)
        for t2_ in tempos[: i + 1]:
            H2 = H(t2_)
            comut = H1 @ H2 - H2 @ H1
            Omega2 += 0.5 * (-1j / hbar)**2 * comut * (dt * dt)
    return Omega1 + Omega2


def evolucao_magnus(
    H: Callable[[float], np.ndarray],
    t0: float,
    t1: float,
    n_passos: int = 64,
    hbar: float = 1.0,
) -> np.ndarray:
    """U(t1,t0) = exp(Ω) com Ω de Magnus ordem 2."""
    Omega = serie_magnus_ordem2(H, t0, t1, n_passos, hbar)
    return expm(Omega)


def fidelidade_unitaria(U: np.ndarray) -> float:
    return float(np.linalg.norm(U.conj().T @ U - np.eye(U.shape[0]), ord="fro"))


def exemplo_campo_magnetico_pulsado(
    B0: float,
    tau: float,
    direcao: str = "x",
) -> Callable[[float], np.ndarray]:
    """
    H(t) = B0 exp(−t²/(2τ²)) · σ  (pulso Gaussiano de campo magnético).
    """
    from .grupos_classicos import geradores_pauli
    sig = geradores_pauli()
    idx = {"x": 0, "y": 1, "z": 2}[direcao]

    def H(t: float) -> np.ndarray:
        return B0 * np.exp(-t**2 / (2 * tau**2)) * sig[idx]

    return H


def comparacao_preservacao_norma(
    H: Callable[[float], np.ndarray],
    psi0: np.ndarray,
    t0: float,
    t1: float,
    n_passos: int = 50,
) -> Dict:
    """
    Compara integrador de Lie-Euler (Magnus-like) com Euler ordinário.
    Ilustra Fig. 12.1 do tratado.
    """
    dt = (t1 - t0) / n_passos
    # Lie
    psi = psi0.copy().astype(complex)
    normas_lie = [1.0]
    for n in range(n_passos):
        t = t0 + n * dt
        U = expm(-1j * H(t) * dt)
        psi = U @ psi
        normas_lie.append(float(np.linalg.norm(psi)))
    # Euler
    psi = psi0.copy().astype(complex)
    normas_euler = [1.0]
    for n in range(n_passos):
        t = t0 + n * dt
        psi = psi - 1j * dt * (H(t) @ psi)
        normas_euler.append(float(np.linalg.norm(psi)))
    return {
        "erro_lie": float(np.max(np.abs(np.array(normas_lie) - 1))),
        "erro_euler": float(np.max(np.abs(np.array(normas_euler) - 1))),
        "referencia": "Wilcke 2026, Cap. 12 §12.4 Fig. 12.1",
    }
