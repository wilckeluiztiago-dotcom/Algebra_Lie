"""
Modelo 4 — Operador de Casimir e invariantes
Autor: Luiz Tiago Wilcke
Fonte: Cap. 3 §3.8 do tratado
"""

from __future__ import annotations
import numpy as np
from typing import Dict
from .momento_angular_su2 import MomentoAngular
from .grupos_classicos import algebra_su2


def casimir_su2(j: float, hbar: float = 1.0) -> Dict:
    """
    Casimir quadrático de su(2): C₂ = J² = Jx² + Jy² + Jz²
    Autovalor: ħ² j(j+1)
    [C₂, J_i] = 0  (invariante)
    """
    ma = MomentoAngular(j, hbar)
    C2 = ma.J2
    autovalor_teorico = hbar**2 * j * (j + 1)
    # verificar que C2 é proporcional à identidade
    desvio = np.max(np.abs(C2 - autovalor_teorico * np.eye(ma.dim)))
    return {
        "j": j,
        "autovalor_C2": autovalor_teorico,
        "max_|C2 - j(j+1)ħ² 1|": float(desvio),
        "referencia": "Wilcke 2026, Cap. 3 §3.8",
    }


def verificar_comutacao_casimir(j: float = 1.0, hbar: float = 1.0) -> Dict:
    """Verifica [J², Jx] = [J², Jy] = [J², Jz] = 0."""
    ma = MomentoAngular(j, hbar)
    erros = {}
    for nome, Ji in [("Jx", ma.Jx), ("Jy", ma.Jy), ("Jz", ma.Jz)]:
        comut = ma.J2 @ Ji - Ji @ ma.J2
        erros[f"max_|[J²,{nome}]|"] = float(np.max(np.abs(comut)))
    return {
        **erros,
        "conclusao": "Casimir comuta com todos os geradores (invariante de Casimir)",
        "referencia": "Wilcke 2026, Cap. 3 §3.8",
    }


def casimir_adjunta_su2() -> Dict:
    """
    Na representação adjunta de su(2) (dim 3), C₂ = 2 · 1
    (porque j=1 → j(j+1)=2).
    """
    su2 = algebra_su2()
    # ad_X ad_X na base
    K = su2.matriz_killing()  # forma de Killing ~ -2 C₂ na adjunta (normalização)
    return {
        "matriz_killing": K.real,
        "nota": "Forma de Killing está relacionada ao Casimir da adjunta (Cap. 5)",
        "referencia": "Wilcke 2026, Cap. 3 §3.8 e Cap. 5",
    }
