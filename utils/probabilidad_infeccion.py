# probabilidad_infección_campo_medio.py
# Autores: José Carlos Riego y Pablo Rodríguez Soria
# Establece y representa la función que describe la probabilidad de infección de un individuo en función del tiempo para el 
# Modelo SIR de Campo Medio. La finalidad de este programa es verificar que dicha probabilidad está normalizada, no excediendo el valor de 1.

import numpy as np
import matplotlib.pyplot as plt
import os

# Definir rutas
cm_dir = os.path.join('data', 'campo_medio')
red_dir = os.path.join('data', 'discreto_red')

## MODELO DE CAMPO MEDIO
# Parámetros del modelo de campo medio
dt = 3.0  # paso de tiempo MC
b = 0.3 # tasa de contagio
infected_cm = np.loadtxt(os.path.join(cm_dir, 'infected_cm_b_0.3_k_0.07_N_10000_dt_3.00.txt'), delimiter=',')  # Vector con la proporción de infectados

P = dt*b*infected_cm    # Probabildiad de infección (proporcional a la proporción de infectados)
# Generamos t basado en la longitud de los datos para evitar errores de dimensión
t = np.arange(len(P)) * dt
tmax_cm = t[-1] 
print(len(P), len(t))
print(infected_cm)

# Gráfica de la evolución temporal de la probabilidad de infección
plt.figure(figsize=(10, 6))
plt.scatter(t, P)
plt.axhline(y=1, color='r', linestyle='-')  # Dibuja una línea horizontal roja en y=1
plt.xlabel('Tiempo (días)')
interval = 10  # Intervalo de 20 días para los ticks en el eje x
plt.xticks(ticks=np.arange(0, tmax_cm + interval, interval))  # Ajustar para que se muestren los ticks cada 10 días
plt.ylabel('Probabilidad de infección')
plt.xlim(0, tmax_cm)  # Establecer el límite del eje x
plt.ylim(0, 1)

# Actualización del título para incluir b y k truncado a dos decimales
title = "Probabilidad de infección en función del tiempo (Campo Medio)"
plt.title(title)

plt.legend(loc='best')
plt.grid()
plt.show()

## MODELO DISCRETO

# Parámetros del modelo discreto
dt = 6.0 # paso de tiempo MC (debe coincidir con nombre del archivo)
b = 0.3 # tasa de contagio
neighbors = 4  # Número de vecinos (cambiar a 8 si se desea tener en cuenta los segundos vecinos)
p_4_neighbors = ((0.22 * np.log(b) + 0.72)) # Probabilidad de contagio en el contacto con un vecinos (4 VECINOS)

# Del ajuste logarítmico de b y p con 4 y 8 vecinos, extraemos la siguiente relación, donde p es la probabilidad de contagio
# en el contacto con un vecino en función del número de vecinos 'neighbors' a tener en cuenta
p = (4 / neighbors) * p_4_neighbors

infected_neighbors_count = np.loadtxt(os.path.join(red_dir, 'infected_neighbors_count_red_b_0.3_k_0.07_N_10000_dt_6.00_neigh_4.txt'), delimiter=',')  # Vector con el número de vecinos infectados en un sitio en medio de la red
P = 1 - (1 - p) ** infected_neighbors_count   # Probabilidad de infección

# Generamos t basado en la longitud de los datos para evitar errores de dimensión
t = np.arange(len(P)) * dt
tmax_red = t[-1]

print(len(P), len(t))
print(infected_neighbors_count)

# Gráfica de la evolución temporal de la probabilidad de infección
plt.figure(figsize=(10, 6))
plt.scatter(t, P)
plt.axhline(y=1, color='r', linestyle='-')  # Dibuja una línea horizontal roja en y=1
plt.xlabel('Tiempo (días)')
interval = 10  # Intervalo de 20 días para los ticks en el eje x
plt.xticks(ticks=np.arange(0, tmax_red + interval, interval))  # Ajustar para que se muestren los ticks cada 10 días
plt.ylabel('Probabilidad de infección')
plt.xlim(0, tmax_red)  # Establecer el límite del eje x
plt.ylim(0, 1)

# Actualización del título para incluir b y k truncado a dos decimales
title = "Probabilidad de infección en función del tiempo (Modelo discreto)"
plt.title(title)

plt.legend(loc='best')
plt.grid()
plt.show()