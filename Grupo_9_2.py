''' 
Universidade Federal do Rio de Janeiro
Escola de Química - EQE 358
Grupo 9
Prazo para entrega: 24/09/2025
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
from numpy.polynomial import Polynomial

# Função para interpolação de Lagrange

def lagrange_interpol(x_vals, y_vals, x_star):
    n = len(x_vals)
    y_star = 0.0
    for i in range(n):
        p = 1.0
        for j in range(n):
            if j != i:
                p *= (x_star - x_vals[j]) / (x_vals[i] - x_vals[j])
        y_star += p * y_vals[i]
    return y_star

# Dados fornecidos

x_data = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
y_data = np.array([0.070225, 0.095651, 0.163071, 0.318481, 0.673870])

# Função real f(x) = sinh(5x)/(x*sinh(5)), mas com lim x ~> 0 f(x) = 5/sinh(5)

def f(x):
    if x == 0:
        return 5.0/np.sinh(5.0)
    else:
        return np.sinh(5*x)/(x*np.sinh(5.0)) 

# O menor DRE do grupo

dres = [123477015, 117205612, 120125499, 122156307, 122090323, 121123945, 120135567]
dre_menor = min(dres)

# Três últimos dígitos do menor DRE do grupo

DRE3 = int(str(dre_menor)[-3:])

# Pontos de teste para interpolação/extrapolação

x_test = [0.0, 0.05, 0.2, 0.5, 0.6, 0.7, 0.8, 0.95, 1.0, DRE3/1000]

# Cálculo dos valores interpolados/extrapolados

results = []
for x_star in x_test:
    y_star = lagrange_interpol(x_data, y_data, x_star)
    results.append((x_star, y_star))

# Exibição dos resultados

print("Resultados da interpolação/extrapolação:")
for x_star, y_star in results:
    print(f"x* = {x_star:.3f}, y* = {y_star:.6f}")

# Polinômio gerado e seus coeficientes
coef = np.polyfit(x_data, y_data, 4)


p = Polynomial(coef[::-1])  
print("\nExpressão do polinômio gerado:")
print(p)

# Gráfico da função real f(x), do polinômio interpolador Pn(x) e do erro Rn(x)

x_plot = np.linspace(0, 1, 101)
y_real = np.array([f(x) for x in x_plot])
y_poly = [lagrange_interpol(x_data, y_data, xx) for xx in x_plot]

erro = y_real - np.array(y_poly)

# Plotagem
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,6))

# Polinômio Pn vs. função real f(x)
ax1.plot(x_data, y_data, 'ro', label='Dados da tabela')
ax1.plot([x for x,_ in results], [y for _,y in results], 'bs', label='Valores interpolados da tabela')
ax1.plot(x_plot, y_real, 'g-', label='Função f(x)')
ax1.plot(x_plot, y_poly, 'b--', label='Polinômio interpolador Pn(x)')

ax1.set_title("Interpolação de Lagrange")
ax1.set_xlabel("x"); ax1.set_ylabel("f(x) e Pn(x)")
ax1.grid(True); ax1.legend()

ax1.text(
    0.025, 0.785, f'Pn(x): {p}', 
    transform=ax1.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7)
)


ax1.text(
    0.025, 0.70, f'f(x) = sinh(5x)/(x·sinh(5))',
    transform=ax1.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7)
)

ax1.text(
    0.025, 0.64, f'f(0) = 5/sinh(5)',
    transform=ax1.transAxes, fontsize=8, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7)
)

# Erro da interpolação Rn(x)
ax2.plot(x_plot, erro, 'm-', label="Função erro Rn(x) = f(x) - Pn(x)")
ax2.axhline(0, color='k', linestyle='--')
ax2.set_title("Erro da Interpolação")
ax2.set_xlabel("x"); ax2.set_ylabel("Rn(x)")
ax2.grid(True); ax2.legend()

# Ajuste de layout e exibição
plt.tight_layout()
plt.show()
