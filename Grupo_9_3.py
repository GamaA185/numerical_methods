''' 
Universidade Federal do Rio de Janeiro
Escola de Química - EQE 358
Grupo 9
Prazo para entrega: 14/10/2025
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

# Bibliotecas necessárias

import numpy as np
import matplotlib.pyplot as plt

# Expressão cujo define o reator CSTR do 3º problema dos Exercícios Computacionais (2025/2)

def f(x):
    return (x + 0.05 - (0.01 * np.exp((15*x)/(1+x))) / (1 + (0.01 * np.exp((15*x)/(1+x)))))

# Função h(a,b) específica do método de Wegstein

''' 
    Função h(a,b) visando manter lambda entre 0 e 1, em virtude do embasamento geométrico 
obtido na similaridade de triângulos no método da secante e no método de Wegstein:

f(a) / (lambda * (b-a)) = -f(b) / (1 - lambda) * (b-a)
=>  lambda = f(a)/(f(a)-f(b))
'''

def h(a, b):
    fa, fb = f(a), f(b)
    return fa / (fa - fb)

# Ajuste de intervalo (só se for necessário, ou seja, f(a) e f(b) tiverem o mesmo sinal ou tu tiver duas raizes muito próximas...)

def encontrar_intervalo_valido(a, b, tentativas=12, passo=0.01):
    fa, fb = f(a), f(b)
    for _ in range(tentativas):
        a -= passo
        b += passo
        fa, fb = f(a), f(b)
        if fa * fb < 0:
            print(f"Intervalo ajustado para: [{a:.3f}, {b:.3f}]")
            return a, b
    print("Nenhum intervalo válido encontrado.")
    return None, None

# Método de Wegstein (Regula Falsi Modificado - Método da Posição Falsa)

def wegstein(a, b, eps=1e-6, delta=1e-8, kmax=100):
    fa = f(a)
    fb = f(b)

    if fa * fb > 0:
        a, b = encontrar_intervalo_valido(a, b)
        if a is None:
            return None, None, 0, None, None
        fa, fb = f(a), f(b)

    k = 0
    while True:

        # lambda ← h(a,b)
        lam = h(a, b)

        # x ← a + lambda * (b - a)
        x = a + lam * (b - a)

        # y ← f(x)

        y = f(x)

        # se y·fa > 0 então
        if y * fa > 0:
            fa = y
            a = x
        else:
            fb = y
            b = x

        # delta ← |b - a|
        delta_x = abs(b - a)
        k += 1

        # condição de parada
        if (delta_x <= eps or abs(y) <= delta) or k >= kmax:
            break

    if k >= kmax:
        print("Número máximo de iterações atingido sem convergência.")
        return None, None, k, a, b

    # retorna raiz
    return x, y, k, a, b

# Aplicação do método para os intervalos estimados incialmente
# Para a fórmula deste reator, não usar x = -1.0, pois f(-1) é indefinido


intervalos = [(-0.09, 0), (0.4, 0.5), (0.8, 0.9)]
resultados = []

for (a, b) in intervalos:
    raiz, fx, it, ai, bi = wegstein(a, b)
    if raiz is not None:
        resultados.append((raiz, fx, it, ai, bi))

# Tabela cujo reporta o valor da raiz, f(x) na raiz [f(raiz)], número de iterações e os valores das duas estimativas iniciais [a,b]

print("\n================= RESULTADOS ================= ")
print("Raiz:\t\tf(x):\t\tIterações:\tIntervalo [a, b]:")
for r in resultados:
    print(f"{r[0]:.6f}\t{r[1]:.2e}\t{r[2]}\t\t[{r[3]:.2f}, {r[4]:.2f}]")

# Cálculo do ponto DRE3/1000

dres = [123477015, 117205612, 120125499, 122156307, 122090323, 121123945, 120135567]
media = int(np.mean(dres))   # média inteira
dre3 = int(str(media)[-3:])  # últimos 3 dígitos
x_dre = dre3 / 1000          # ponto DRE3/1000

# Gráfico com todas as raízes e f(DRE3/1000)

x_vals = np.linspace(-0.9999, 1.0000, 600)
y_vals = f(x_vals)

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(x_vals, y_vals, label='f(x)', color='blue')
ax.axhline(0, color='black', lw=0.8)

for r in resultados:
    ax.plot(r[0], 0, 'ro', label=f'Raiz = {r[0]:.3f}')

ax.plot(x_dre, f(x_dre), 'ks', label=f'f({x_dre:.3f}) = {f(x_dre):.2e}')
ax.set_title('Gráfico da função f(x) do reator CSTR via método de Wegstein')
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.legend()
ax.grid(True)

# Cria a tabela com os resultados

tabela_dados = [[f"{r[0]:.6f}", f"{r[1]:.2e}", r[2], f"[{r[3]:.2f}, {r[4]:.2f}]"] for r in resultados]
colunas = ["Raiz", "f(x)", "Iterações", "Intervalo estimado"]

# Adiciona a tabela abaixo do gráfico

ax.table(cellText=tabela_dados, colLabels=colunas, cellLoc='center', loc='bottom', bbox=[0, -0.7, 1, 0.3])
plt.subplots_adjust(left=0.1, bottom=0.45)
fig.text(0.5, 0.32,
         "Tabela 1 – Resultados numéricos obtidos pelo método de Wegstein para o reator CSTR.",
         ha='center', va='center', fontsize=9, style='italic')

plt.show()
