# SIR_determinista.py
# Autores: José Carlos Riego y Pablo Rodríguez Soria
# Modelo SIR determinista que obtiene la evolución temporal de la propagación de un virus.

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os

# Definir rutas
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), 'data')
det_dir = os.path.join(data_dir, 'determinista')
os.makedirs(det_dir, exist_ok=True)

# Parámetros del modelo
b = 0.3     # tasa de transmisión
k = 1/14    # tasa de recuperación
tmax = 100  # tiempo máximo de integración de las EDOs

# Condiciones iniciales
I0 = 0.01  # pequeña fracción de infectados iniciales
R0 = 0.0  # inicialmente no hay recuperados
S0 = 1 - I0 - R0  # el resto de la población es susceptible

# Función que define las EDOs del modelo SIR
def deriv(y, t, b, k):
    S, I, R = y # noqa: E741
    dSdt = -b * S * I
    dIdt = b * S * I - k * I
    dRdt = k * I
    return dSdt, dIdt, dRdt

# Vector de tiempo (días)
t = np.linspace(0, tmax, tmax)

# Condiciones iniciales para la integración
y0 = S0, I0, R0

# Integración de las EDOs
ret = odeint(deriv, y0, t, args=(b, k))
S, I, R = ret.T  # noqa: E741

# Guardamos los vectores con el conteo de infectados, sanos y recuperados
# Nota: En el determinista N está normalizado a 1 (proporciones), por lo que no lo incluimos en el nombre
base_name = f'determinista_b_{b}_k_{k:.2f}'
np.savetxt(os.path.join(det_dir, f'susceptible_{base_name}.txt'), S, fmt='%f', delimiter=',')
np.savetxt(os.path.join(det_dir, f'infected_{base_name}.txt'), I, fmt='%f', delimiter=',')
np.savetxt(os.path.join(det_dir, f'recovered_{base_name}.txt'), R, fmt='%f', delimiter=',')
np.savetxt(os.path.join(det_dir, f'time_{base_name}.txt'), t, fmt='%f', delimiter=',')

# Gráfica de resultados
plt.figure(figsize=(10, 6))
plt.plot(t, S, 'b', label='Susceptibles')
plt.plot(t, I, 'r', label='Infectados')
plt.plot(t, R, 'g', label='Recuperados')
plt.xlabel('Tiempo (días)')
interval = 10  # Intervalo de 20 días para los ticks en el eje x
plt.xticks(ticks=np.arange(0, tmax + interval, interval))  # Ajustar para que se muestren los ticks cada 10 días
plt.ylabel('Proporción de la población')

# Actualización del título para incluir b y k truncado a dos decimales
title = f"Modelo SIR determinista (b = {b}, 1/k = {(1/k):.2f})"
plt.title(title)

plt.legend(loc='best')
plt.grid()
plt.show()

