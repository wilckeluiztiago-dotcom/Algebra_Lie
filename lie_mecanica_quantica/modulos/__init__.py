"""
Pacote de Álgebras de Lie aplicadas à Mecânica Quântica
Autor: Luiz Tiago Wilcke
Fonte: Tratado "Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais" (2026)
"""

from .algebra_lie_basica import AlgebraDeLie, verificar_identidade_jacobi
from .grupos_classicos import algebra_su2, algebra_heisenberg, geradores_pauli
from .representacoes import coeficientes_clebsch_gordan_su2

from .oscilador_heisenberg import OsciladorHarmonico, espectro_oscilador
from .momento_angular_su2 import MomentoAngular, autovalores_j2_jz
from .adicao_spins_cg import adicionar_spins, base_acoplada
from .casimir_invariantes import casimir_su2, verificar_comutacao_casimir
from .efeito_zeeman import HamiltonianoZeeman, espectro_zeeman
from .sistema_dois_niveis import HamiltonianoRabi, frequencia_rabi
from .estados_coerentes import estado_coerente, valor_esperado_coerente
from .evolucao_magnus_mq import evolucao_magnus, serie_magnus_ordem2
from .stone_von_neumann import representacao_schrodinger, verificar_relacao_ccr
from .hidrogenio_angular import HarmonicosEsfericos, energia_angular

__version__ = "1.0.0"
__autor__ = "Luiz Tiago Wilcke"
