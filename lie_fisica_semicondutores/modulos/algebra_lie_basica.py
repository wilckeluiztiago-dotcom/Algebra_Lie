"""
Módulo: algebra_lie_basica
Implementa a definição axiomática, colchete de Lie, identidade de Jacobi,
constantes de estrutura, representação adjunta e forma de Killing.

Referências no tratado de Luiz Tiago Wilcke (2026):
- Capítulo 1: Definição Axiomática das Álgebras de Lie
- §1.3 Definição Axiomática e Propriedades Fundamentais
- §1.5 O Problema de Estrutura e as Constantes de Estrutura
- §1.7 Interpretações Profundas da Identidade de Jacobi
- §1.8 O Homomorfismo e Representação Unificada Adjunta
- Capítulo 5: Critério de Cartan e a Forma de Killing
"""

from __future__ import annotations
import numpy as np
from typing import List, Tuple, Optional, Callable


class AlgebraDeLie:
    """
    Álgebra de Lie finito-dimensional sobre K = R ou C.

    Axiomas (Definição 1.3 do tratado):
    1. Espaço vetorial de dimensão finita.
    2. Colchete bilinear [ , ] : g × g → g
    3. Antissimetria: [X,Y] = -[Y,X]
    4. Identidade de Jacobi: [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0
    """

    def __init__(
        self,
        dimensao: int,
        constantes_estrutura: np.ndarray,
        nome: str = "g",
        base_nomes: Optional[List[str]] = None,
    ):
        """
        Parâmetros
        ----------
        dimensao : int
            dim g
        constantes_estrutura : array (n,n,n)
            c^k_{ij} tais que [e_i, e_j] = sum_k c^k_{ij} e_k
            (convenção do tratado §1.5)
        nome : str
            nome simbólico da álgebra
        base_nomes : list[str]
            nomes dos geradores da base
        """
        self.dimensao = dimensao
        self.constantes_estrutura = np.asarray(constantes_estrutura, dtype=complex)
        self.nome = nome
        self.base_nomes = base_nomes or [f"e{i}" for i in range(dimensao)]

        if self.constantes_estrutura.shape != (dimensao, dimensao, dimensao):
            raise ValueError(
                f"constantes_estrutura deve ter forma ({dimensao},{dimensao},{dimensao})"
            )

        # Verificação automática dos axiomas
        self._verificar_antissimetria()
        if not verificar_identidade_jacobi(self.constantes_estrutura):
            raise ValueError(
                f"As constantes de estrutura de {nome} violam a identidade de Jacobi "
                "(tratado §1.7)."
            )

    def colchete_lie(self, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        """
        Calcula [X,Y] = sum_{i,j,k} X^i Y^j c^k_{ij} e_k
        (definição via constantes de estrutura – Cap. 1 §1.5)
        """
        X = np.asarray(X, dtype=complex).ravel()
        Y = np.asarray(Y, dtype=complex).ravel()
        if X.shape[0] != self.dimensao or Y.shape[0] != self.dimensao:
            raise ValueError("Vetores devem ter dimensão igual a dim g")

        resultado = np.zeros(self.dimensao, dtype=complex)
        for i in range(self.dimensao):
            for j in range(self.dimensao):
                for k in range(self.dimensao):
                    resultado[k] += (
                        X[i] * Y[j] * self.constantes_estrutura[i, j, k]
                    )
        return resultado

    def _verificar_antissimetria(self) -> None:
        """Verifica c^k_{ij} = -c^k_{ji} (axioma de antissimetria)."""
        for i in range(self.dimensao):
            for j in range(self.dimensao):
                if not np.allclose(
                    self.constantes_estrutura[i, j, :],
                    -self.constantes_estrutura[j, i, :],
                ):
                    raise ValueError(
                        f"Antissimetria violada em índices ({i},{j}) – Cap. 1 §1.3"
                    )

    def matriz_adjunta(self, X: np.ndarray) -> np.ndarray:
        """
        Representação adjunta ad_X : Y ↦ [X,Y]
        (tratado §1.8 – Homomorfismo e Representação Unificada Adjunta)
        """
        X = np.asarray(X, dtype=complex).ravel()
        matriz = np.zeros((self.dimensao, self.dimensao), dtype=complex)
        for j in range(self.dimensao):
            ej = np.zeros(self.dimensao, dtype=complex)
            ej[j] = 1.0
            matriz[:, j] = self.colchete_lie(X, ej)
        return matriz

    def forma_de_killing(self, X: np.ndarray, Y: np.ndarray) -> complex:
        """
        Forma de Killing κ(X,Y) = Tr(ad_X ∘ ad_Y)
        (Capítulo 5 – Critério de Cartan e a Forma de Killing)
        """
        adX = self.matriz_adjunta(X)
        adY = self.matriz_adjunta(Y)
        return np.trace(adX @ adY)

    def matriz_killing(self) -> np.ndarray:
        """Matriz da forma de Killing na base escolhida."""
        K = np.zeros((self.dimensao, self.dimensao), dtype=complex)
        for i in range(self.dimensao):
            ei = np.zeros(self.dimensao)
            ei[i] = 1.0
            for j in range(self.dimensao):
                ej = np.zeros(self.dimensao)
                ej[j] = 1.0
                K[i, j] = self.forma_de_killing(ei, ej)
        return K

    def eh_semissimples(self, tolerancia: float = 1e-10) -> bool:
        """
        Critério de Cartan: g é semissimples ⇔ det(κ) ≠ 0
        (tratado Cap. 5 §5.4)
        """
        det = np.linalg.det(self.matriz_killing())
        return abs(det) > tolerancia

    def __repr__(self) -> str:
        return f"AlgebraDeLie({self.nome}, dim={self.dimensao})"


def verificar_identidade_jacobi(
    constantes: np.ndarray, tolerancia: float = 1e-10
) -> bool:
    """
    Verifica a identidade de Jacobi para todas as ternas de geradores.
    Interpretação homológica (Cap. 1 §1.7.2): d² = 0 no complexo de Chevalley-Eilenberg.
    """
    n = constantes.shape[0]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                # [ei, [ej, ek]] + cic
                termo = np.zeros(n, dtype=complex)
                # [ej, ek]
                comut_jk = constantes[j, k, :]
                # [ei, [ej,ek]]
                for m in range(n):
                    termo += comut_jk[m] * constantes[i, m, :]
                # ciclando
                comut_ki = constantes[k, i, :]
                for m in range(n):
                    termo += comut_ki[m] * constantes[j, m, :]
                comut_ij = constantes[i, j, :]
                for m in range(n):
                    termo += comut_ij[m] * constantes[k, m, :]

                if np.any(np.abs(termo) > tolerancia):
                    return False
    return True


def forma_de_killing(algebra: AlgebraDeLie, X: np.ndarray, Y: np.ndarray) -> complex:
    """Função auxiliar – mesma da classe."""
    return algebra.forma_de_killing(X, Y)


def representacao_adjunta(algebra: AlgebraDeLie, X: np.ndarray) -> np.ndarray:
    """Função auxiliar – matriz de ad_X."""
    return algebra.matriz_adjunta(X)
