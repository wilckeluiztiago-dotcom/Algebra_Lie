# Álgebras de Lie aplicadas à Física de Semicondutores

**Autor do projeto e do tratado de referência:** Luiz Tiago Wilcke  
**Tratado fonte:** *Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais* (1ª edição, 2026)  
**Departamento de Matemática e Física Teórica**

Este repositório implementa, em Python + NumPy/SciPy, aplicações completas e executáveis de álgebras de Lie a problemas centrais da física de semicondutores. Todas as construções matemáticas estão explicitamente referenciadas aos capítulos e seções do tratado anexo.

## Origem das ideias no livro

| Conceito / Modelo | Capítulo / Seção do tratado | Aplicação em semicondutores |
|-------------------|-----------------------------|-----------------------------|
| Colchete de Lie, Jacobi, constantes de estrutura | Cap. 1 §1.3–1.5, §1.7 | Estrutura algébrica dos geradores de simetria efetivos |
| Álgebra de Heisenberg \(\mathfrak{h}_3\) | Cap. 1 §1.6.3 e Cap. 10 §10.2 | Operadores de posição/momento na aproximação de massa efetiva e níveis de Landau |
| \(\mathfrak{su}(2)\), momento angular, Clebsch-Gordan | Cap. 1 §1.6.1–1.6.2, Cap. 3 §3.7, Cap. 10 §10.3–10.4 | Spin, acoplamento spin-órbita (Rashba/Dresselhaus), estrutura fina de éxcitons |
| Representação adjunta e forma de Killing | Cap. 1 §1.8 e Cap. 5 | Invariantes e classificação de termos de quebra de simetria |
| Sistemas de raízes e diagramas de Dynkin | Cap. 6–7 | Classificação de simetrias contínuas efetivas |
| Fórmula de Weyl e pesos | Cap. 9 | Multiplicidades e degenerescências de bandas |
| Série de Magnus e integradores em grupos de Lie | Cap. 12 §12.2–12.4 | Evolução temporal unitária sob campos elétricos/magnéticos dependentes do tempo |
| Matrix Product States e simetrias de Lie | Cap. 11 | Métodos de tensores para heteroestruturas e pontos quânticos correlacionados |

## Equações fundamentais utilizadas

### 1. Álgebra de Heisenberg e massa efetiva

Na aproximação de massa efetiva o Hamiltoniano de envelope é

$$
H = \frac{\mathbf{p}^2}{2m^*} + V(\mathbf{r}),
$$

onde os operadores de posição e momento satisfazem a álgebra de Heisenberg (tratado Cap. 1 §1.6.3 e Cap. 10 §10.2)

$$
[x_i,p_j] = i\hbar\delta_{ij},\qquad [x_i,x_j]=[p_i,p_j]=0.
$$

### 2. Níveis de Landau

Na presença de um campo magnético \(\mathbf{B}=B\hat{z}\) o momento canônico é substituído pelo momento cinético \(\boldsymbol{\pi}=\mathbf{p}+e\mathbf{A}\). Os operadores

$$
a = \frac{\ell_B}{\sqrt{2}\hbar}(\pi_x - i\pi_y),\qquad
a^\dagger = \frac{\ell_B}{\sqrt{2}\hbar}(\pi_x + i\pi_y)
$$

obedecem \([a,a^\dagger]=1\) (álgebra de Heisenberg). Os níveis de Landau são

$$
E_n = \hbar\omega_c\Bigl(n+\frac12\Bigr),\qquad \omega_c=\frac{eB}{m^*},\qquad n=0,1,2,\dots
$$

### 3. Acoplamento spin-órbita de Rashba

O termo de Rashba em um poço quântico 2DEG (interface assimétrica) é

$$
H_R = \alpha_R\bigl(\sigma_x k_y - \sigma_y k_x\bigr) = \alpha_R\boldsymbol{\sigma}\cdot(\mathbf{k}\times\hat{z}),
$$

onde \(\boldsymbol{\sigma}\) são as matrizes de Pauli (geradores de \(\mathfrak{su}(2)\), tratado Cap. 1 §1.6.1 e Cap. 10 §10.3). O espectro é

$$
E_\pm(\mathbf{k}) = \frac{\hbar^2k^2}{2m^*}\pm\alpha_R k.
$$

### 4. Acoplamento de Dresselhaus

Em cristais sem centro de inversão (GaAs, InAs, etc.) aparece o termo cúbico de Dresselhaus, que na aproximação linear (poços estreitos) reduz-se a

$$
H_D = \beta_D\bigl(\sigma_x k_x - \sigma_y k_y\bigr).
$$

A combinação Rashba + Dresselhaus é descrita por um Hamiltoniano de \(\mathfrak{su}(2)\) efetivo

$$
H_{\mathrm{SO}} = \boldsymbol{\Omega}(\mathbf{k})\cdot\boldsymbol{\sigma},
$$

com vetor de precessão \(\boldsymbol{\Omega}(\mathbf{k})\) linear em \(\mathbf{k}\).

### 5. Evolução temporal – série de Magnus

Para um Hamiltoniano dependente do tempo \(H(t)\) (campo elétrico ou magnético pulsado) a evolução unitária é

$$
U(t) = \mathcal{T}\exp\Bigl(-\frac{i}{\hbar}\int_0^t H(s)\,ds\Bigr) = \exp\bigl(\Omega(t)\bigr),
$$

onde a série de Magnus (tratado Cap. 12 §12.2) lê-se

$$
\Omega(t) = \Omega_1 + \Omega_2 + \Omega_3 + \cdots,
$$

$$
\Omega_1 = -\frac{i}{\hbar}\int_0^t H(t_1)\,dt_1,
$$

$$
\Omega_2 = \frac12\Bigl(-\frac{i}{\hbar}\Bigr)^2\int_0^t dt_1\int_0^{t_1}dt_2\,[H(t_1),H(t_2)],
$$

e assim por diante. O truncamento preserva unitariedade (integrador geométrico em grupo de Lie).

### 6. Estrutura fina de éxcitons – Clebsch-Gordan

O elétron (\(j_e\)) e o buraco (\(j_h\)) acoplam-se segundo

$$
\mathbf{J} = \mathbf{j}_e + \mathbf{j}_h.
$$

Os estados \(\lvert J,M\rangle\) são obtidos pelos coeficientes de Clebsch-Gordan (tratado Cap. 3 §3.7 e Cap. 10 §10.4)

$$
\lvert J,M\rangle = \sum_{m_e+m_h=M}\langle j_e m_e;j_h m_h\lvert J M\rangle\,\lvert j_e m_e\rangle\otimes\lvert j_h m_h\rangle.
$$

### 7. Graus de liberdade de vale (grafeno / TMD)

Em grafeno ou dicalcogenetos de metais de transição o Hamiltoniano efetivo de dois vales (K e K′) é isomorfo a um \(\mathfrak{su}(2)\) de vale:

$$
H_{\mathrm{vale}} = v_F\boldsymbol{\sigma}\cdot\mathbf{k} + \tau_z\Delta\,\sigma_z,
$$

onde \(\tau_z=\pm 1\) rotula o vale. Em sistemas de quatro componentes (spin × vale) aparece um \(\mathfrak{su}(4)\) efetivo.

## Estrutura do projeto

```
lie_fisica_semicondutores/
├── README.md
├── requisitos.txt
├── modulos/
│   ├── __init__.py
│   ├── algebra_lie_basica.py      # Cap. 1 e 5
│   ├── grupos_classicos.py        # su(2), su(3), Heisenberg
│   ├── representacoes.py          # Clebsch-Gordan, Weyl
│   ├── sistemas_raizes.py         # raízes e Dynkin
│   ├── niveis_landau.py           # Heisenberg + campo B
│   ├── spin_orbita.py             # Rashba, Dresselhaus, su(2)
│   ├── evolucao_magnus.py         # série de Magnus (Cap. 12)
│   ├── excitons_cg.py             # estrutura fina via CG
│   └── vales_su2.py               # Hamiltoniano de vale
├── exemplos/
│   ├── demo_landau.py
│   ├── demo_rashba_dresselhaus.py
│   ├── demo_magnus_pulso.py
│   └── demo_exciton.py
└── testes/
    └── test_semicondutores.py
```

## Instalação e execução

```bash
cd lie_fisica_semicondutores
pip install -r requisitos.txt
python -m testes.test_semicondutores
python -m exemplos.demo_landau
python -m exemplos.demo_rashba_dresselhaus
python -m exemplos.demo_magnus_pulso
python -m exemplos.demo_exciton
```

## Filosofia (Prefácio do tratado)

> “A filosofia deste livro é que intuição geométrica e rigor algébrico são complementares, não opostos.”  
> — Luiz Tiago Wilcke, Março de 2026

Cada módulo começa pela motivação física do semicondutor e só depois formaliza os objetos de álgebra de Lie, exatamente como recomendado no tratado.

## Citação

```
Wilcke, Luiz Tiago. Álgebras de Lie — Simetria, Estrutura e Métodos Computacionais.
Tratado de Pós-Graduação em Matemática e Física Teórica, 1ª ed., 2026.
```

---
*Autor: Luiz Tiago Wilcke.*
