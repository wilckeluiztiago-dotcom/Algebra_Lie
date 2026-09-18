"""
Pacote de Álgebras de Lie aplicadas à Física de Partículas
Autor: Luiz Tiago Wilcke
Fonte: Tratado "Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais" (2026)
"""

from .algebra_lie_basica import (
    AlgebraDeLie,
    verificar_identidade_jacobi,
    forma_de_killing,
    representacao_adjunta,
)
from .grupos_classicos import (
    algebra_su2,
    algebra_su3,
    algebra_so3,
    algebra_heisenberg,
    geradores_pauli,
    geradores_gell_mann,
)
from .sistemas_raizes import (
    SistemaDeRaizes,
    matriz_de_cartan,
    diagrama_dynkin_simples,
)
from .representacoes import (
    coeficientes_clebsch_gordan_su2,
    pesos_representacao_su3,
    formula_dimensao_weyl,
)
from .fisica_particulas import (
    modelo_padrao,
    decomposicao_mesons,
    decomposicao_barions,
    grupo_gut_su5,
    cadeia_excepcional_e8,
)

__version__ = "1.0.0"
__autor__ = "Luiz Tiago Wilcke"
