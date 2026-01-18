# probabilidad_recuperación.py
# Autores: José Carlos Riego y Pablo Rodríguez Soria
# Establece y representa la función que describe la probabilidad de recuperación de un individuo frente al virus en función del tiempo.

import numpy as np
import matplotlib.pyplot as plt

# Valor del tiempo típico tau
k = 1/14  # Tasa de recuperación (población recuperada al 50%)
tmax = 100  # tiempo máximo de integración de las EDOs
tau = (1 / (k * np.log(2))) # Relación entre k y tau. Tau es el tiempo t para el cual (1 - 1/e) de la población recuperada

# Función para generar la probabilidad de recuperación P(x)
def recovery_function(x,tc):
    return (1 - np.exp(-x/tc))

# Creamos una gama de valores de x de 0 a 30 días (1 mes)
x_values = np.arange(0, tmax)

# Calculamos P(x) para el rango de x usando la función que encontramos
P_x = recovery_function(x_values, tau)

# Gráfico de P(x)
plt.figure(figsize=(10, 6))
plt.plot(x_values, P_x, label="Probabilidad de Recuperación P(x)")
title = f"Probabilidad de recuperación (tau = {tau:.2f},  1/k = {(1/k):.2f})"
plt.title(title)
plt.xlabel("Tiempo desde la infección (días)")
plt.ylabel("Probabilidad de Recuperación")
plt.xlim(0, tmax)
plt.ylim(0, 1.05)
plt.legend()
plt.grid(True)
plt.show()
