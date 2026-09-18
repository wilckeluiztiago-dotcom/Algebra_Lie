"""
Módulo: evolucao_magnus
Série de Magnus e integradores geométricos em grupos de Lie
para evolução temporal em semicondutores sob campos dependentes do tempo.

Referências no tratado de Luiz Tiago Wilcke (2026):
- Capítulo 12: Sistemas Dinâmicos e Simulações
  §12.2 A Série de Magnus
  §12.4 Integradores Simpléticos em Grupos de Lie
- Cap. 2 §2.8  Fórmula de Baker-Campbell-Hausdorff
"""

from __future__ import annotations
import numpy as np
from typing import Callable, Tuple, List, Optional
from scipy.linalg import expm


def serie_magnus_ordem2(
    H: Callable[[float], np.ndarray],
    t0: float,
    t1: float,
    n_passos: int = 64,
) -> np.ndarray:
    """
    Aproximação de Magnus até segunda ordem (tratado Cap. 12 §12.2):

        Ω₁ = −(i/ħ) ∫_{t0}^{t1} H(t) dt
        Ω₂ = (1/2) (−i/ħ)² ∫∫_{t0≤t2≤t1} [H(t1), H(t2)] dt1 dt2

    (aqui trabalhamos em unidades ħ = 1; o fator −i é absorvido na definição
    de H hermitiano → gerador anti-hermitiano).

    Retorna o gerador Ω ≈ Ω₁ + Ω₂ (matriz anti-hermitiana).
    """
    dt = (t1 - t0) / n_passos
    tempos = np.linspace(t0 + dt / 2, t1 - dt / 2, n_passos)

    # Ω₁
    Omega1 = np.zeros_like(H(t0), dtype=complex)
    for t in tempos:
        Omega1 += -1j * H(t) * dt

    # Ω₂ (aproximação de ponto médio)
    Omega2 = np.zeros_like(Omega1)
    for i, t1_ in enumerate(tempos):
        H1 = H(t1_)
        for t2_ in tempos[: i + 1]:
            H2 = H(t2_)
            comut = H1 @ H2 - H2 @ H1
            Omega2 += (-1j)**2 * 0.5 * comut * (dt * dt)

    return Omega1 + Omega2


def evolucao_unitaria_magnus(
    H: Callable[[float], np.ndarray],
    t0: float,
    t1: float,
    n_passos: int = 64,
) -> np.ndarray:
    """
    Operador de evolução U(t1,t0) ≈ exp(Ω) com Ω da série de Magnus
    truncada na ordem 2. Preserva unitariedade por construção
    (integrador geométrico – tratado Cap. 12 §12.4).
    """
    Omega = serie_magnus_ordem2(H, t0, t1, n_passos)
    return expm(Omega)


def integrador_lie_euler(
    H: Callable[[float], np.ndarray],
    psi0: np.ndarray,
    t0: float,
    t1: float,
    n_passos: int = 100,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrador de Lie-Euler (ordem 1) em grupo de Lie:

        ψ_{n+1} = exp(−i H(t_n) Δt) ψ_n

    Mantém a norma ‖ψ‖ = 1 a menos de erro de arredondamento.
    (Comparar com a Figura 12.1 do tratado: integrador geométrico vs Euler ordinário.)
    """
    dt = (t1 - t0) / n_passos
    tempos = np.linspace(t0, t1, n_passos + 1)
    psi = psi0.astype(complex).copy()
    trajetoria = [psi.copy()]

    for n in range(n_passos):
        t = tempos[n]
        U_passo = expm(-1j * H(t) * dt)
        psi = U_passo @ psi
        trajetoria.append(psi.copy())

    return tempos, np.array(trajetoria)


def exemplo_pulso_eletrico(
    E0: float,
    tau: float,
    m_efetiva: float,
    alpha_R: float = 0.0,
    dim: int = 2,
) -> Callable[[float], np.ndarray]:
    """
    Hamiltoniano dependente do tempo de um elétron 2DEG sob um pulso
    elétrico Gaussiano + termo de Rashba constante:

        H(t) = p²/(2m*) + e E(t)·x  +  α_R (σ×k)·ẑ

    Para ilustração usamos uma matriz 2×2 efetiva no espaço de spin
    com um termo de Zeeman pulsado:

        H(t) = (Δ0 + E0 exp(−t²/(2τ²))) σ_z  +  α_R σ_x
    """
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    def H(t: float) -> np.ndarray:
        envelope = E0 * np.exp(-t**2 / (2 * tau**2))
        return (envelope) * sigma_z + alpha_R * sigma_x

    return H


def fidelidade_unitaria(U: np.ndarray) -> float:
    """
    Medida de quão unitária é a matriz: ‖U†U − 1‖_F.
    Deve ser ~0 para o integrador de Magnus/Lie.
    """
    return float(np.linalg.norm(U.conj().T @ U - np.eye(U.shape[0]), ord="fro"))


def comparacao_euler_vs_magnus(
    H: Callable[[float], np.ndarray],
    psi0: np.ndarray,
    t0: float,
    t1: float,
    n_passos: int = 50,
) -> dict:
    """
    Compara a preservação da norma entre Euler ordinário (que sai da
    variedade) e o integrador de Lie-Euler (que permanece no grupo).
    Ilustra a Figura 12.1 do tratado.
    """
    # Lie-Euler
    tempos, traj_lie = integrador_lie_euler(H, psi0, t0, t1, n_passos)
    normas_lie = np.array([np.linalg.norm(psi) for psi in traj_lie])

    # Euler ordinário (diferença finita na equação de Schrödinger)
    dt = (t1 - t0) / n_passos
    psi = psi0.astype(complex).copy()
    normas_euler = [np.linalg.norm(psi)]
    for n in range(n_passos):
        t = t0 + n * dt
        psi = psi - 1j * dt * (H(t) @ psi)   # Euler explícito
        normas_euler.append(np.linalg.norm(psi))
    normas_euler = np.array(normas_euler)

    return {
        "tempos": tempos,
        "normas_lie": normas_lie,
        "normas_euler": normas_euler,
        "erro_norma_lie": float(np.max(np.abs(normas_lie - 1.0))),
        "erro_norma_euler": float(np.max(np.abs(normas_euler - 1.0))),
        "referencia": "Wilcke 2026, Cap. 12 §12.4 e Fig. 12.1",
    }
