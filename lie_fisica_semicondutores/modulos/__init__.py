"""
Pacote de Álgebras de Lie aplicadas à Física de Semicondutores
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
    algebra_heisenberg,
    geradores_pauli,
)
from .representacoes import (
    coeficientes_clebsch_gordan_su2,
    formula_dimensao_weyl,
)
from .niveis_landau import (
    HamiltonianosLandau,
    espectro_landau,
    operadores_ladder_landau,
)
from .spin_orbita import (
    HamiltonianosSpinOrbita,
    espectro_rashba,
    espectro_rashba_dresselhaus,
)
from .evolucao_magnus import (
    serie_magnus_ordem2,
    evolucao_unitaria_magnus,
    integrador_lie_euler,
)
from .excitons_cg import (
    estados_exciton_cg,
    hamiltoniano_estrutura_fina,
)
from .vales_su2 import (
    HamiltonianosVale,
    espectro_grafeno_massivo,
)

__version__ = "1.0.0"
__autor__ = "Luiz Tiago Wilcke"
