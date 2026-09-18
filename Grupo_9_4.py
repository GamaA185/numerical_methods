''' 
Universidade Federal do Rio de Janeiro
Escola de Química - EQE 358
Grupo 9
Prazo para entrega: 02/11/2025
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
import pandas as pd
import matplotlib.pyplot as plt

# Cálculo de Phi a partir dos três últimos dígitos do maior DRE do grupo

DREs = [123477015, 117205612, 120125499, 122156307, 122090323, 121123945, 120135567]
DRE = max(DREs)
DRE3 = int(str(DRE)[-3:])
Phi = DRE3 / 200.0  # Phi = DRE3 / 200

# Dados de entrada da função a ser criada

a, b = 0.0, 1.0
N = 2
delta = 1e-6
hmin = 1e-8
eps = hmin

# Função integranda 

def f(x, Phi=Phi):
    return 3 * x * (np.sinh(Phi * x) / np.sinh(Phi))

# Implementação da Regra de Simpson Composta com Extrapolação de Richardson

def simpson_composta_richardson_fiel(f, a, b, N, delta, hmin, Phi):

    # ETAPA 1 – Cálculo inicial

    S0 = f(a, Phi) + f(b, Phi)
    h = (b - a) / N

    js_impar = np.arange(1, N, 2)
    Simpar = np.sum(f(a + js_impar * h, Phi)) if js_impar.size > 0 else 0.0

    if N > 1:
        js_par = np.arange(2, N, 2)
        Spar = np.sum(f(a + js_par * h, Phi)) if js_par.size > 0 else 0.0
    else:
        Spar = 0.0

    I = (h / 3.0) * (S0 + 4 * Simpar + 2 * Spar)
    history = [(1, N, h, I)]

    # ETAPA 2 – Processo Recursivo

    while True:
        Ivelho = I
        N += N       
        h *= 0.5       

        Spar = Spar + Simpar  

        js_novos = np.arange(1, N, 2)
        Simpar = np.sum(f(a + js_novos * h, Phi))

        I = (h / 3.0) * (S0 + 4 * Simpar + 2 * Spar)
        history.append((2, N, h, I))

        if abs(I - Ivelho) <= delta:
            convergiu = True
            break
        if abs(h) <= eps:
            convergiu = False
            break

    # ETAPA 3 – Extrapolação de Richardson

    IR = (16 * I - Ivelho) / 15.0
    history.append((3, N, h, IR))

    return I, IR, N, history, convergiu

# Execução

I, IR, N_final, history, convergiu = simpson_composta_richardson_fiel(
    f, a, b, N, delta, hmin, Phi
)

# Criação e exibição da tabela com resultados finais
# Criar DataFrame apenas com I, IR e N

df_resultados = pd.DataFrame({
    "Descrição": ["Integral (antes de Richardson)", "Integral (após Richardson)", "Número de intervalos usados"],
    "Valor": [I, IR, N_final]
})

# Criação e exibição das tabelas:

# Criar DataFrame com resultados finais 

df_resultados = pd.DataFrame({
    "Descrição": [
        "Integral (I)", 
        "Integral (IR)", 
        "Número de intervalos (N)"
    ],
    "Valor": [
        I,
        IR,
        N_final
    ]
})

# Criar DataFrame para o histórico de cálculos 

df_hist = pd.DataFrame({
    "Etapa": [h[0] for h in history],
    "Intervalos (N)": [h[1] for h in history],
    "h": [h[2] for h in history],
    "Integral": [h[3] for h in history]
})

# Preparar df_resultados formatado

df_resultados_display = df_resultados.copy()
def fmt_val(v):
    if isinstance(v, (float, np.floating)):
        return f"{v:.12e}"
    return str(v)
df_resultados_display['Valor'] = df_resultados_display['Valor'].apply(fmt_val)

# Preparar histórico formatado

df_hist_display = df_hist.copy()
df_hist_display['h'] = df_hist_display['h'].map(lambda x: f"{x:.3e}")
df_hist_display['Integral'] = df_hist_display['Integral'].map(lambda x: f"{x:.12e}")

# Criar figura com duas tabelas

fig, axes = plt.subplots(2, 1, figsize=(10, 8))
for ax in axes:
    ax.axis('off')

# Tabela resultados finais

table_res = axes[0].table(
    cellText=df_resultados_display.values,
    colLabels=df_resultados_display.columns,
    cellLoc='center',
    loc='center'
)
table_res.auto_set_font_size(False)
table_res.set_fontsize(10)
table_res.scale(1, 1.4)
axes[0].set_title("Resultados Finais", pad=10)

# Tabela histórico

table_hist = axes[1].table(
    cellText=df_hist_display.values,
    colLabels=df_hist_display.columns,
    cellLoc='center',
    loc='center'
)
table_hist.auto_set_font_size(False)
table_hist.set_fontsize(9)
table_hist.scale(1, 1.2)
axes[1].set_title("Histórico de iterações – Regra 1/3 de Simpson Composta", pad=10)

plt.tight_layout()
plt.show()

# Impressão final no console

print("\n===== RESULTADOS FINAIS =====")
print(f"I (antes de Richardson) = {I:.12e}")
print(f"IR (após Richardson)   = {IR:.12e}")
print(f"2N (número de intervalos) = {N_final}")


# ================================================== Q.E.D. ================================================== #
