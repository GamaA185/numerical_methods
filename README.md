# Projetos de Métodos Numéricos - EQE 358

Este repositorio reune cinco scripts em Python desenvolvidos para exercicios computacionais da disciplina EQE 358, na Escola de Quimica da UFRJ. Os programas aplicam metodos numéricos em problemas de aproximação de funções, interpolação, busca de raízes, integração numérica e solução de sistemas de equações diferenciais ordinárias.

## Requisitos

- Python 3.10 ou superior
- NumPy
- Matplotlib
- Pandas

Instalação das dependencias:

```bash
pip install numpy matplotlib pandas
```

## Como executar

Cada arquivo pode ser executado diretamente:

```bash
python GRUPO_9_1.py
python Grupo_9_2.py
python Grupo_9_3.py
python Grupo_9_4.py
python Grupo_9_5.py
```

Os scripts imprimem resultados numéricos no terminal e, em geral, exibem gráficos com Matplotlib.

## Arquivos

### `GRUPO_9_1.py`

Aproxima a funcao `exp(-x^2)` por serie de Taylor com criterio de convergencia baseado na precisao da maquina. O programa:

- calcula o epsilon da maquina;
- define `delta = sqrt(epsilon)`;
- avalia a aproximacao em pontos fixos e em um ponto baseado no maior DRE do grupo;
- compara os valores obtidos com `numpy.exp`;
- calcula erros absoluto e relativo;
- gera um grafico comparando a funcao exata e a aproximacao.

### `Grupo_9_2.py`

Implementa interpolacao polinomial pelo metodo de Lagrange para dados tabelados. O programa:

- calcula valores interpolados e extrapolados em pontos de teste;
- monta o polinomio interpolador de grau 4;
- compara o polinomio com a funcao analitica `sinh(5x)/(x*sinh(5))`;
- trata o limite da funcao em `x = 0`;
- plota o polinomio, a funcao real e a funcao erro.

### `Grupo_9_3.py`

Resolve uma equacao nao linear associada a um modelo de reator CSTR usando uma forma do metodo de Wegstein, equivalente a uma posicao falsa modificada. O programa:

- define a funcao nao linear do problema;
- procura raizes em tres intervalos iniciais;
- ajusta intervalos quando nao ha troca de sinal;
- imprime raiz, valor de `f(x)`, numero de iteracoes e intervalo final;
- gera grafico da funcao com as raizes encontradas.

### `Grupo_9_4.py`

Calcula uma integral definida usando a regra 1/3 de Simpson composta com refinamento sucessivo e extrapolacao de Richardson. O programa:

- define o parametro `Phi` a partir do maior DRE do grupo;
- integra a funcao `3*x*sinh(Phi*x)/sinh(Phi)` no intervalo `[0, 1]`;
- refina o numero de subintervalos ate atingir a tolerancia definida;
- aplica extrapolacao de Richardson;
- exibe tabelas com resultados finais e historico de iteracoes.

### `Grupo_9_5.py`

Simula um sistema dinamico de biomassa e substrato em um reator usando Adams-Bashforth de terceira ordem, inicializado com Runge-Kutta de terceira ordem. O programa:

- define parametros do modelo biologico;
- calcula a taxa especifica de crescimento `mu(x2)`;
- resolve o sistema para diferentes valores de vazao volumetrica;
- plota a evolucao temporal da biomassa;
- plota a evolucao temporal do substrato.
