''' 
Universidade Federal do Rio de Janeiro
Escola de Química - EQE 358
Grupo_9_1.py

Docente responsável: D.Sc. Prof. Argimiro R. Secchi (PEQ/COPPE/UFRJ)

Discentes do grupo 9:
    Andrew de Carvalho Leite Gama - 123477015
    Bruno Lyra Ferreira - 117205612
    Geórgia de Moura Veras - 120125499
    Giulia de Andrade Lanciano - 122156307
    Janine Aparecida Pereira Gazoni - 122090323
    Pedro Lucas Caetano de Oliveira - 121123945
    Thiago Ferreira Campos - 120135567
    '''

import numpy as np
import matplotlib.pyplot as plt

def epsilon():
    """
    Calcula a precisão da máquina.

    Retorna:
        A precisão da máquina (epsilon).
    """
    eps = 1.0
    while eps + 1.0 != 1.0:
        eps *= 0.5
    eps *= 2.0
    return eps

def exp2T(x, delta):
    """
    Aproxima a função f(x) = exp(-x^2) usando a série de Taylor, com
    critério de convergência delta.

    Argumentos:
        x (float): O valor de x para o qual a função será aproximada.
        delta (float): O critério de convergência.

    Retorna:
        Uma tupla contendo o valor aproximado da função, o grau do polinômio (n) e o valor de y.
    """
    y = x**2
    m = 0

    # Redução da magnitude de y para melhorar a convergência
    while y >= 1:
        m += 1
        y /= 2

    n = 0
    T = 1.0
    S = 1.0

    # Cálculo da série de Taylor
    while True:
        n += 1
        T = -y * T / n
        S += T
        if abs(T / S) <= delta:
            break

    # Restauração da magnitude de S
    for _ in range(m):
        S *= S

    return S, n, y

def main():
    """
    Função principal que executa os testes, calcula os erros
    e gera o gráfico comparativo.
    """
    # Lista de DREs 
    dres = [123477015, 120125499, 120135567, 117205612, 122156307, 121123945, 122090323]

    # Seleção do maior DRE
    maior_dre = max(dres)

    # Pegar os três últimos dígitos como string, mantendo zeros à esquerda
    ultimos_tres_str = str(maior_dre)[-3:]

    # Converter para inteiro (mantém valor correto mesmo com zeros à esquerda)
    ultimos_tres_int = int(ultimos_tres_str)

    # Calcular x_dre
    x_dre = ultimos_tres_int / 100

    # Valores de x para teste
    test_x = [0.1, 1.0, 2.0, x_dre]

    #Precisão da máquina, eps, e critério de convergência, raiz quadrada de eps

    eps = epsilon()
    delta = np.sqrt(eps)

    print(f"Precisão da máquina (epsilon): {eps:.16f}")
    print(f"Critério de convergência (delta): {delta:.16f}")
    print("-" * 50)

    for x in test_x:
        # Calcular a aproximação
        f_approx, n, _ = exp2T(x, delta)

        # Calcular o valor exato usando a biblioteca numpy
        f_exact = np.exp(-x**2)

        # Calcular os erros
        erro_absoluto = abs(f_exact - f_approx)
        erro_relativo = erro_absoluto / abs(f_exact) if f_exact != 0 else 0

        print(f"Resultados para x = {x}:")
        print(f"  Valor aproximado: {f_approx:.16f}")
        print(f"  Valor exato (Numpy): {f_exact:.16f}")
        print(f"  Grau do polinômio (n): {n}")
        print(f"  Erro absoluto: {erro_absoluto:.16f}")
        print(f"  Erro relativo: {erro_relativo:.16f}")
        print("-" * 50)

    # --- Geração do gráfico ---
    # Gerar pontos para a curva da função exata
    x_vals = np.linspace(-2, 2, 400)
    y_exact = np.exp(-x_vals**2)

    # Gerar pontos para a curva da função aproximada
    y_approx = [exp2T(val, delta)[0] for val in x_vals]

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_exact, label='Função Exata $f(x) = e^{-x^2}$ (Numpy)', color='blue')
    plt.plot(x_vals, y_approx, label='Função Aproximada (Taylor)', color='red', linestyle='--')
    plt.title('Comparação entre a Função Exata e a Aproximação de Taylor')
    plt.xlabel('x')
    plt.ylabel('$f(x)$')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
