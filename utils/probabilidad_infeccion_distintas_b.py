# probabilidad_infección_disntitas_b.py
# Autores: José Carlos Riego y Pablo Rodríguez Soria
# Representa la probabilidad de infección de un individuo en función del tiempo en el Modelo SIR de Campo Medio para distintos valores de b. 
# La finalidad de este programa es verificar que dicha probabilidad está normalizada, no excediendo el valor de 1, para distintos valores de b.

import numpy as np
import matplotlib.pyplot as plt
import os

# Definir rutas
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), 'data', 'other', 'calculo_p(b)')


## MODELO DE CAMPO MEDIO
# Parámetros del modelo de campo medio
dt = 3  # paso de tiempo MC
b = 0.3 # tasa de contagio

infected_cm_0_1 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.1.txt'), delimiter=',')  # Vector con la proporción de infectados

P_0_1 = dt*b*infected_cm_0_1    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_2 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.2.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_2 = dt*b*infected_cm_0_2    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_3 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.3.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_3 = dt*b*infected_cm_0_3    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_4 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.4.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_4 = dt*b*infected_cm_0_4    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_5 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.5.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_5 = dt*b*infected_cm_0_5    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_6 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.6.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_6 = dt*b*infected_cm_0_6    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_7 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.7.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_7 = dt*b*infected_cm_0_7    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_8 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.8.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_8 = dt*b*infected_cm_0_8    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_9 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_0.9.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_9 = dt*b*infected_cm_0_9    # Probabildiad de infección (proporcional a la proporción de infectados)

infected_cm_0_11 = np.loadtxt(os.path.join(data_dir, 'infected_cm_b_1.0.txt'), delimiter=',')  # Vector con la proporción de infectados
P_0_11 = dt*b*infected_cm_0_11    # Probabildiad de infección (proporcional a la proporción de infectados)

tmax = 100  # tiempo máximo de cálculo de las proporciones SIR

# Generar t basado en la longitud real de los datos para evitar errores de dimensión
t = np.arange(len(P_0_1)) * dt

# Gráfica de la evolución temporal de la probabilidad de infección
plt.figure(figsize=(10, 6))
plt.plot(t, P_0_1, 'o:', label = 'b = 0.1')
plt.plot(t, P_0_2, 'o:', label = 'b = 0.2')
plt.plot(t, P_0_3, 'o:', label = 'b = 0.3')
plt.plot(t, P_0_4, 'o:', label = 'b = 0.4')
plt.plot(t, P_0_5, 'o:', label = 'b = 0.5')
plt.plot(t, P_0_6, 'o:', label = 'b = 0.6')
plt.plot(t, P_0_7, 'o:', label = 'b = 0.7')
plt.plot(t, P_0_8, 'o:', label = 'b = 0.8')
plt.plot(t, P_0_9, 'o:', label = 'b = 0.9')
plt.plot(t, P_0_11, 'o:', label = 'b = 1.0')
plt.axhline(y=1, color='r', linestyle='-')  # Dibuja una línea horizontal roja en y=1
plt.xlabel('Tiempo (días)')
interval = 10  # Intervalo de 20 días para los ticks en el eje x
plt.xticks(ticks=np.arange(0, tmax + interval, interval))  # Ajustar para que se muestren los ticks cada 10 días
plt.ylabel('Probabilidad de infección')
plt.xlim(0, tmax)  # Establecer el límite del eje x
plt.ylim(0, 1.1)

# Actualización del título para incluir b y k truncado a dos decimales
title = "Probabilidad de infección en función del tiempo (Campo Medio)"
plt.title(title)

plt.legend(loc='best')

plt.show()

## MODELO DISCRETO

# Parámetros del modelo discreto
dt = 5 # paso de tiempo MC
b = 0.3 # tasa de contagio
neighbors = 8  # Número de vecinos (cambiar a 8 si se desea tener en cuenta los segundos vecinos)
p_4_neighbors = ((0.22 * np.log(b) + 0.72)) # Probabilidad de contagio en el contacto con un vecinos (4 VECINOS)

# Del ajuste logarítmico de b y p con 4 y 8 vecinos, extraemos la siguiente relación, donde p es la probabilidad de contagio
# en el contacto con un vecino en función del número de vecinos 'neighbors' a tener en cuenta
p = (4 / neighbors) * p_4_neighbors

infected_neighbors_count = np.loadtxt(os.path.join(data_dir, 'infected_neighbors_count.txt'), delimiter=',')  # Vector con el número de vecinos infectados en un sitio en medio de la red
P = 1 - (1 - p) ** infected_neighbors_count   # Probabilidad de infección

# Generar t basado en la longitud real de los datos para evitar errores de dimensión
t = np.arange(len(P)) * dt
tmax_red = t[-1]

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