# Álgebras de Lie aplicadas à Física de Partículas

**Autor do projeto e do tratado de referência:** Luiz Tiago Wilcke  
**Tratado fonte:** *Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais* (1ª edição, 2026)  
**Departamento de Matemática e Física Teórica**

Este repositório implementa, em Python puro + NumPy, vários modelos concretos de física de partículas a partir dos conceitos desenvolvidos no tratado anexo. Todas as ideias, notações e teoremas utilizados estão explicitamente referenciados aos capítulos e seções do livro.

## Origem das ideias no livro (referências precisas)

| Conceito implementado | Capítulo / Seção do livro | Observação |
|-----------------------|---------------------------|------------|
| Definição axiomática, colchete, Jacobi, constantes de estrutura | Cap. 1 (Def. 1.1–1.3, §1.5, §1.7) | Base de toda a biblioteca |
| Álgebra de Heisenberg \(\mathfrak{h}_3\) e oscilador harmônico | Cap. 1 §1.6.3 e Cap. 10 §10.2 | Princípio de incerteza e quantização |
| \(\mathfrak{su}(2)\), momento angular, Clebsch-Gordan | Cap. 1 §1.6.1–1.6.2, Cap. 3 §3.7, Cap. 10 §10.3–10.4 | Composição de spins |
| Forma de Killing e semissimplicidade | Cap. 5 | Critério de Cartan |
| Sistemas de raízes, matriz de Cartan, diagramas de Dynkin | Cap. 6 e 7 | Classificação |
| Representações de \(\mathfrak{su}(3)\), pesos, fórmula de Weyl | Cap. 9 e Cap. 10 §10.6.1 | Octeto e decupleto de hádrons |
| Teoria de gauge Yang-Mills | Cap. 10 §10.5 | Tensor de curvatura \(F_{\mu\nu}\) |
| **Modelo Padrão** \(SU(3)_c\times SU(2)_L\times U(1)_Y\) | Cap. 10 §10.6 | Dimensões dos geradores → 8+3+1 bósons |
| Mésons \(3\otimes\bar{3}=1\oplus 8\) e bárions \(3\otimes 3\otimes 3\) | Cap. 10 §10.6.1 | Predição do \(\Omega^-\) |
| GUTs: \(SU(5)\), \(SO(10)\), cadeia \(E_8\supset\dots\supset\) SM | Cap. 8 §8.5–8.6 e Cap. 10 §10.7 | Dimensões e decaimento do próton |
| Excepcionais \(G_2,F_4,E_6,E_7,E_8\) | Cap. 8 | Dados estruturais e Dynkin |

Todas as variáveis de código estão em **português** (colchete_lie, constantes_estrutura, geradores, etc.), seguindo a nomenclatura do tratado.

## Estrutura do projeto

```
lie_fisica_particulas/
├── README.md
├── requisitos.txt
├── modulos/
│   ├── __init__.py
│   ├── algebra_lie_basica.py      # Cap. 1 – colchete, Jacobi, adjunta, Killing
│   ├── grupos_classicos.py        # su(2), su(3), so(3), heisenberg
│   ├── sistemas_raizes.py         # Cap. 6–7 – raízes, Cartan, Dynkin
│   ├── representacoes.py          # Cap. 3 e 9 – pesos, Clebsch-Gordan, Weyl
│   └── fisica_particulas.py       # Cap. 10 – SM, Yang-Mills, GUTs, hádrons
├── exemplos/
│   ├── demo_momento_angular.py
│   ├── demo_modelo_padrao.py
│   ├── demo_sabor_su3.py
│   └── demo_guts.py
└── testes/
    └── test_basicos.py
```

## Instalação rápida

```bash
cd lie_fisica_particulas
pip install -r requisitos.txt
python -m exemplos.demo_modelo_padrao
```

## Filosofia (extraída do Prefácio do tratado)

> “A filosofia deste livro é que intuição geométrica e rigor algébrico são complementares, não opostos.”  
> — Luiz Tiago Wilcke, Março de 2026

Os módulos começam sempre com a motivação geométrica/física (fluxos, gauge, spins) e só depois formalizam os objetos algébricos, exatamente como o tratado recomenda.

## Licença e citação

Se utilizar este código ou o tratado, cite:

```
Wilcke, Luiz Tiago. Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais.
Tratado de Pós-Graduação em Matemática e Física Teórica, 1ª ed., 2026.
```

---
*Gerado a partir da análise completa do PDF anexo (220 páginas).*
