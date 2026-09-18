# Álgebras de Lie aplicadas à Mecânica Quântica

**Autor do projeto e do tratado de referência:** Luiz Tiago Wilcke  
**Tratado fonte:** *Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais* (1ª edição, 2026)  
**Departamento de Matemática e Física Teórica**

Este repositório implementa **10 modelos completos e executáveis** de Mecânica Quântica construídos a partir das álgebras de Lie desenvolvidas no tratado. Cada modelo extrai explicitamente definições, teoremas e interpretações do livro (principalmente Capítulos 1, 3, 10 e 12).

## Origem das ideias no livro

| # | Modelo | Capítulo / Seção do tratado |
|---|--------|-----------------------------|
| 1 | Oscilador harmônico (álgebra de Heisenberg) | Cap. 1 §1.6.3 e Cap. 10 §10.2 |
| 2 | Momento angular e \(\mathfrak{su}(2)\) | Cap. 1 §1.6.1–1.6.2 e Cap. 10 §10.3 |
| 3 | Adição de spins (Clebsch-Gordan) | Cap. 3 §3.7 e Cap. 10 §10.4 |
| 4 | Operador de Casimir e invariantes | Cap. 3 §3.8 |
| 5 | Efeito Zeeman (spin em campo magnético) | Cap. 10 §10.3 |
| 6 | Sistema de dois níveis / precessão de Rabi | Cap. 1 §1.6 e Cap. 10 |
| 7 | Estados coerentes de Glauber | Cap. 10 §10.2 (Heisenberg) |
| 8 | Evolução temporal — série de Magnus | Cap. 12 §12.2–12.4 |
| 9 | Teorema de Stone–von Neumann | Cap. 10 §10.1 |
| 10 | Átomo de hidrogênio (parte angular \(\mathfrak{so}(3)\)) | Cap. 1 §1.6.2 e Cap. 10 §10.3 |

## Equações fundamentais

### 1. Álgebra de Heisenberg e oscilador harmônico

$$
[x,p] = i\hbar,\qquad a = \sqrt{\frac{m\omega}{2\hbar}}x + \frac{i}{\sqrt{2m\omega\hbar}}p,\qquad [a,a^\dagger]=1
$$

$$
H = \hbar\omega\Bigl(a^\dagger a + \tfrac12\Bigr),\qquad E_n = \hbar\omega\Bigl(n+\tfrac12\Bigr)
$$

(Tratado Cap. 1 §1.6.3 e Cap. 10 §10.2)

### 2. Momento angular — \(\mathfrak{su}(2)\)

$$
[J_i,J_j] = i\hbar\varepsilon_{ijk}J_k,\qquad \mathbf{J}^2 = \hbar^2 j(j+1),\qquad J_z = \hbar m
$$

(Tratado Cap. 1 §1.6.1–1.6.2 e Cap. 10 §10.3)

### 3. Clebsch-Gordan

$$
|j_1 m_1\rangle\otimes|j_2 m_2\rangle = \sum_{J,M}\langle j_1 m_1;j_2 m_2|JM\rangle\,|JM\rangle
$$

$$
V^{j_1}\otimes V^{j_2} = \bigoplus_{J=|j_1-j_2|}^{j_1+j_2} V^J
$$

(Tratado Cap. 3 §3.7 e Cap. 10 §10.4, eq. 10.13)

### 4. Operador de Casimir

$$
C_2 = \sum_a T^a T^a \quad\text{(quadrático)},\qquad [C_2,T^b]=0
$$

Para \(\mathfrak{su}(2)\): \(C_2 = \mathbf{J}^2\).

(Tratado Cap. 3 §3.8)

### 5. Efeito Zeeman

$$
H_Z = -\boldsymbol{\mu}\cdot\mathbf{B} = \frac{g\mu_B}{\hbar}\mathbf{J}\cdot\mathbf{B}
$$

(Tratado Cap. 10 §10.3)

### 6. Hamiltoniano de Rabi (dois níveis)

$$
H = \frac{\hbar\omega_0}{2}\sigma_z + \frac{\hbar\Omega}{2}(\sigma_+ e^{-i\omega t} + \sigma_- e^{i\omega t})
$$

(Tratado Cap. 1 e Cap. 10 — estrutura \(\mathfrak{su}(2)\))

### 7. Estados coerentes de Glauber

$$
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty\frac{\alpha^n}{\sqrt{n!}}|n\rangle,\qquad a|\alpha\rangle = \alpha|\alpha\rangle
$$

(Tratado Cap. 10 §10.2)

### 8. Série de Magnus

$$
U(t) = \exp\bigl(\Omega(t)\bigr),\qquad
\Omega_1 = -\frac{i}{\hbar}\int_0^t H(t_1)\,dt_1,\qquad
\Omega_2 = \frac12\Bigl(-\frac{i}{\hbar}\Bigr)^2\int_0^t dt_1\int_0^{t_1}[H(t_1),H(t_2)]\,dt_2
$$

(Tratado Cap. 12 §12.2)

### 9. Stone–von Neumann

Toda representação irredutível da álgebra de Heisenberg em espaço de Hilbert (com regularidade) é unitariamente equivalente à representação de Schrödinger.

(Tratado Cap. 10 §10.1)

### 10. Parte angular do hidrogênio

$$
H_{\mathrm{ang}} = \frac{\mathbf{L}^2}{2mr^2},\qquad Y_{lm}(\theta,\phi)\ \text{autovetores de }\mathbf{L}^2\text{ e }L_z
$$

(Tratado Cap. 1 §1.6.2 e Cap. 10 §10.3)

## Estrutura do projeto

```
lie_mecanica_quantica/
├── README.md
├── requisitos.txt
├── modulos/
│   ├── __init__.py
│   ├── algebra_lie_basica.py
│   ├── grupos_classicos.py
│   ├── representacoes.py
│   ├── oscilador_heisenberg.py      # Modelo 1
│   ├── momento_angular_su2.py       # Modelo 2
│   ├── adicao_spins_cg.py           # Modelo 3
│   ├── casimir_invariantes.py       # Modelo 4
│   ├── efeito_zeeman.py             # Modelo 5
│   ├── sistema_dois_niveis.py       # Modelo 6
│   ├── estados_coerentes.py         # Modelo 7
│   ├── evolucao_magnus_mq.py        # Modelo 8
│   ├── stone_von_neumann.py         # Modelo 9
│   └── hidrogenio_angular.py        # Modelo 10
├── exemplos/
│   ├── demo_oscilador.py
│   ├── demo_spins_cg.py
│   ├── demo_zeeman_rabi.py
│   └── demo_magnus_coerentes.py
└── testes/
    └── test_mq.py
```

## Instalação e execução

```bash
cd lie_mecanica_quantica
pip install -r requisitos.txt
python -m testes.test_mq
python -m exemplos.demo_oscilador
python -m exemplos.demo_spins_cg
python -m exemplos.demo_zeeman_rabi
python -m exemplos.demo_magnus_coerentes
```

## Filosofia (Prefácio do tratado)

> “A filosofia deste livro é que intuição geométrica e rigor algébrico são complementares, não opostos.”  
> — Luiz Tiago Wilcke, Março de 2026

## Citação

```
Wilcke, Luiz Tiago. Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais.
Tratado de Pós-Graduação em Matemática e Física Teórica, 1ª ed., 2026.
```

---
*Gerado a partir da análise completa do PDF anexo (220 páginas).*
