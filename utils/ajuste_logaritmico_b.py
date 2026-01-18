# ajuste_logaritmico_b.py
# Autores: José Carlos Riego y Pablo Rodríguez Soria
# Representa p frente a b, cuyos valores han sido obtenidos mediante la comparación visual entre el modelo continuo y el modelo discreto 
# (en este caso para 8 vecinos), buscando la mayor similitud posible entre las gráficas de las proporciones SIR de ambos modelos,
# y obtiene la función logarítmica que mejor los ajusta, con el fin de obtener p(b). 

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Parámetros
dt = 5  # paso de tiempo de la simulación MC del modelo discreto
b_values = np.array([0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]) # valores de b
p_values = np.array([0.3 , 0.36, 0.45, 0.51, 0.55, 0.59, 0.63, 0.65, 0.68])

# Definimos una función logarítmica para el ajuste
def logarithmic_model(b, a, c):
    return a * np.log(b) + c

# Ajustar la función logarítmica a los datos
popt_log, pcov_log = curve_fit(logarithmic_model, b_values, p_values)

# Parámetros optimizados para el modelo logarítmico
a_opt, c_opt = popt_log

# Crear puntos b para la curva ajustada logarítmica
b_fit_log = np.linspace(b_values.min(), b_values.max(), 100)

# Calcular los valores ajustados de p para el modelo logarítmico
p_fit_log = logarithmic_model(b_fit_log, a_opt, c_opt)

# Graficar los datos originales y la curva ajustada logarítmica
plt.figure(figsize=(10, 6))
plt.scatter(b_values, p_values, color='red', label='Datos originales')
plt.plot(b_fit_log, p_fit_log, 'blue', label='Ajuste Logarítmico')
plt.title(f'b frente a la probabilidad de contagio p y su ajuste logarítmico (dt = {dt})')
plt.xlabel('b')
plt.ylabel('p')
plt.legend()
plt.grid(True)
plt.show()
