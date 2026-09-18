''' 
Universidade Federal do Rio de Janeiro
Escola de Química - EQE 358
Grupo 9
Prazo para entrega: 19/11/2025
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

# Parâmetros

t_final = 15.0       # h
V = 10.0             # m³
x2f = 4.0            # kg/m³
mu0 = 0.53           # h⁻¹
ki = 0.4545          # m³/kg
km = 0.12            # kg/m³
Y = 0.4              # rendimento
h = 0.1              # passo

times = np.arange(0, t_final + 1e-12, h)

# Condições iniciais

x1_0 = 1.0
x2_0 = 0.5

# Fator em comum 

def mu(x2):
    return mu0 * x2 / (km + x2 + ki * x2**2)

# RHS do sistema

def rhs(x, D):
    x1, x2 = x
    m = mu(x2)
    dx1 = m * x1 - D * x1
    dx2 = -(1/Y) * m * x1 - D * x2 + D * x2f
    return np.array([dx1, dx2])

# Runge Kutta de 3º ordem

def rk3_step(x, D, h):
    k1 = rhs(x, D)
    k2 = rhs(x + 0.5 * h * k1, D)
    k3 = rhs(x + h * (-k1 + 2 * k2), D)
    return x + h * ( (1/6) * k1 + (2/3) * k2 + (1/6) * k3 )

# Integração AB3 (com 2 passos iniciais RK3)

def integrate_ab3(D, x0, tvec, h):
    n = len(tvec)
    sol = np.zeros((n, 2))
    sol[0] = x0

    # 1º passo
    sol[1] = rk3_step(sol[0], D, h)
    # 2º passo
    sol[2] = rk3_step(sol[1], D, h)

    # Adams–Bashforth 3º ordem
    for i in range(2, n-1):
        f_n   = rhs(sol[i], D)
        f_nm1 = rhs(sol[i-1], D)
        f_nm2 = rhs(sol[i-2], D)

        sol[i+1] = sol[i] + (h/12)*(23*f_n - 16*f_nm1 + 5*f_nm2)

        # evitar valores negativos numéricos
        sol[i+1] = np.maximum(sol[i+1], 0)

    return sol

# Cenários de fluxo volumétrico F (m³/h)

DREs = [123477015, 117205612, 120125499, 122156307, 122090323, 121123945, 120135567]
DRE = np.random.choice(DREs)
DRE3 = int(str(DRE)[-3:])

scenarios = {
    "F = 3": 3.0,
    "F = 5": 5.0,
    f"F = DRE3/300 ={float(DRE3/300):.4f}": DRE3/300
}

# Simulação

results = {}
for name, F in scenarios.items():
    D = F / V
    sol = integrate_ab3(D, np.array([x1_0, x2_0]), times, h)
    results[name] = sol

# Gráfico c/ Biomassa

plt.figure(figsize=(8,5))
for name in results:
    plt.plot(times, results[name][:,0], label=name)

plt.xlabel("Tempo (h)")
plt.ylabel("Biomassa (kg/m³)")
plt.title("Evolução da biomassa")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Gráfico c/ Substrato

plt.figure(figsize=(8,5))
for name in results:
    plt.plot(times, results[name][:,1], label=name)

plt.xlabel("Tempo (h)")
plt.ylabel("Substrato (kg/m³)")
plt.title("Evolução do substrato")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
