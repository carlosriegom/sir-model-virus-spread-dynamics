# SIR_discreto_red.py
# Autores: José Carlos Riego y Pablo Rodríguez Soria
# Modelo SIR discreto de red 2D, que obtiene la evolución temporal de la propagación de un virus, para los casos de 4 y 8 vecinos.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
from scipy.integrate import odeint
import os

# Definir rutas
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), 'data')
videos_dir = os.path.join(os.path.dirname(script_dir), 'videos')
red_dir = os.path.join(data_dir, 'discreto_red')
os.makedirs(red_dir, exist_ok=True)
os.makedirs(videos_dir, exist_ok=True)


# Parámetros
N = 10000  # Asegúrate de que N sea un cuadrado perfecto
side_length = int(np.sqrt(N))  # Determina el tamaño de un lado de la red cuadrada
b = 0.3  # Tasa de contagio
k = 1/14  # Tasa de recuperación (población recuperada al 50%)
tau = (1 / (k * np.log(2))) # Relación entre k y tau. Tau es el tiempo t para el cual (1 - 1/e) de la población recuperada
neighbors = 4  # Número de vecinos (cambiar a 8 si se desea tener en cuenta los segundos vecinos)
print('Número de vecinos por sitio en la red: ', neighbors)
p_4_neighbors = ((0.20 * np.log(b) + 0.70)) # Probabilidad de contagio en el contacto con un vecinos (4 VECINOS)
print(p_4_neighbors)

# Del ajuste logarítmico de b y p con 4 y 8 vecinos, extraemos la siguiente relación, donde p es la probabilidad de contagio
# en el contacto con un vecino en función del número de vecinos 'neighbors' a tener en cuenta
p = (4 / neighbors) * p_4_neighbors  

I0 = 0.01  # Proporción inicial de individuos infectados
S0 = 1 - I0 # Proporción inicial de individuos sanos
side_length = int(np.sqrt(N))  # Determina el tamaño de un lado de la red cuadrada
tmax = 100  # tiempo máximo en días para la simulación
dt = 6 # Incremento de tiempo en días => ESTE ES EL MEJOR dt POSIBLE
steps = int(tmax / dt) + 1 # Número total de pasos de simulación para 160 días

# Contadores para susceptibles, infectados y recuperados
susceptible_count = np.zeros(steps)
infected_count = np.zeros(steps)
recovered_count = np.zeros(steps)

# Inicialización de la lista que contiene el número de vecinos infectados (para el sitio (sqrt(N)/2, sqrt(N)/2 ) en el medio de la red)
infected_neighbors_count = []
infected_neighbors_count = np.zeros(int(tmax / dt))

# Inicialización de la matriz de la red
network = np.ones((side_length, side_length), dtype=int)
initial_infected = np.random.choice(N, int(N * I0), replace=False)
network.ravel()[initial_infected] = 2
infection_time = -np.ones((side_length, side_length), dtype=int)
infection_time.ravel()[initial_infected] = 0

# Función de recuperación
def recovery_function(x,tc):
    return (1 - np.exp(-x/tc))

fist_run = False

# Función para actualizar el estado de la red en cada paso de Monte Carlo
def update_network(frame_num, img, network, infection_time, ax, susceptible_count, infected_count, recovered_count):
    new_network = network.copy()
    # Restablecer contadores para este frame
    S = I = R = 0 # noqa: E741

    global fist_run

    # En el frame_num = 0 (primer frame) tenemos la red inicial
    if frame_num == 0:
        if not fist_run:
            fist_run = True
            infected_count[frame_num] = I0
            susceptible_count[frame_num] = 1 - I0
            recovered_count[frame_num] = 0

            # Contamos los vecinos infectados del sitio en el centro de la red
            i = int(np.sqrt(N) / 2)
            j = int(np.sqrt(N) / 2)
            infected_neighbors = 0

            # Segundos vecinos: horizontal, vertical y diagonal
            second_neighbors = [(-1, -1), (-1, 1), (1, -1), (1, 1)]

            # Revisar los segundos vecinos con condiciones de contorno periódicas
            for dx, dy in second_neighbors:
                nx, ny = (i + dx) % side_length, (j + dy) % side_length
                if new_network[nx, ny] == 2:
                    infected_neighbors += 1

            infected_neighbors_count[frame_num] = infected_neighbors
        
        else:
            # Omitir procesamiento en la segunda llamada al frame 0
            return (img,)

    else: 
    # Proceso de infección y recuperación
        for i in range(side_length):
            for j in range(side_length):
                # Bucle para si el individuo es susceptible de contagio
                if new_network[i, j] == 1:  # Susceptible

                    # Inicializamos el número de vecinos infectados como 0
                    infected_neighbors = 0

                    # Bucle para 4 u 8 vecinos
                    if neighbors == 4 or neighbors == 8:
                        # Vecinos directos
                        direct_neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1)]

                        # Revisar los vecinos directos con condiciones de contorno periódicas
                        for dx, dy in direct_neighbors:
                            nx, ny = (i + dx) % side_length, (j + dy) % side_length
                            if new_network[nx, ny] == 2:
                                infected_neighbors += 1

                        # Bucle para 8 vecinos
                        if neighbors == 8:
                            # Segundos vecinos: horizontal, vertical y diagonal
                            second_neighbors = [(-1, -1), (-1, 1), (1, -1), (1, 1)]

                            # Revisar los segundos vecinos con condiciones de contorno periódicas
                            for dx, dy in second_neighbors:
                                nx, ny = (i + dx) % side_length, (j + dy) % side_length
                                if new_network[nx, ny] == 2:
                                    infected_neighbors += 1

                            # Se coloca en un punto en medio de la red para guardar en una lista el número de vecinos infectados a cada frame
                            if i == np.sqrt(N) / 2 and j == i == np.sqrt(N) / 2: # Punto en medio de la red
                                #infected_neighbors_count.append(infected_neighbors)
                                infected_neighbors_count[frame_num] = infected_neighbors
                                day = frame_num*dt
                                #print(frame_num, day, infected_neighbors_count[frame_num], infected_neighbors)

                    # Probabilidad de contagio
                    infection_prob = 1 - (1 - p) ** infected_neighbors

                    # Posible contagio
                    if np.random.random() < infection_prob:
                        new_network[i, j] = 2  # Infección
                        infection_time[i, j] = frame_num

                # Bucle para si el individuo está infectado
                elif new_network[i, j] == 2:  # Infectado
                    infection_duration = (frame_num - infection_time[i, j]) * dt
                    recovery_prob = recovery_function(infection_duration, tau)

                    # Posible recuperación
                    if np.random.random() < recovery_prob:
                        new_network[i, j] = 3  # Recuperación

        # Contar el estado después de aplicar todos los cambios
        for i in range(side_length):
            for j in range(side_length):
                if new_network[i, j] == 1:
                    S += 1
                elif new_network[i, j] == 2:
                    I += 1 # noqa: E741
                elif new_network[i, j] == 3:
                    R += 1

        # Actualización de la red y de la imagen de la animación
        network[:, :] = new_network
        img.set_array(network)

        # Almacenamiento de las proporciones de cada estado
        susceptible_count[frame_num] = (S / N)
        infected_count[frame_num] = (I / N)
        recovered_count[frame_num] = (R / N)

    # Actualización del título con el día y la semana
    day = frame_num*dt   # Convertir frame_num en días => el frame con el índice 0 es el frame_num = 1 (primer frame)
    week = int(day / 7)  # Calcular la semana actual
    ax.set_title(f'Día {day:.1f} | Semana {week + 1} | N = {N} individuos')
    print(frame_num, day)
    return (img,)

# Crear la figura para la animación y ajustar los ejes
fig, ax = plt.subplots()
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(vmin=1, vmax=3)
img = ax.imshow(network, interpolation='none', norm=norm, cmap=cmap)

# Iniciar la animación
ani = animation.FuncAnimation(fig, update_network, fargs=(img, network, infection_time, ax, susceptible_count, infected_count, recovered_count),
                              frames=steps, interval=200, repeat=False)

# Guardar la animación como un vídeo => descomentar sólo si se va a guardar el vídeo porque da error el resto
#ani.save(os.path.join(videos_dir, f'simulacion_red_b_{b}_k_{k:.2f}_N_{N}_dt_{dt:.2f}_neigh_{neighbors}.mp4'), writer='ffmpeg', fps=3)  # Guardar como video MP4

# Leyenda para el código de colores
colors = [cmap(norm(1)), cmap(norm(2)), cmap(norm(3))]
labels = ['Sano', 'Infectado', 'Recuperado']
patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(3)]
plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.show()

# Después de completar la animación, ajustamos la gráfica para mostrar los datos hasta 160 días
plt.figure(figsize=(10, 6))
time_vector = np.arange(0,tmax+1,dt)
#print(len(time_vector))
print(infected_count)
print(time_vector)

# Si queremos graficar solamente los puntos
plt.scatter(time_vector, susceptible_count, label='Susceptibles', color='blue')
plt.scatter(time_vector, infected_count, label='Infectados', color='red')
plt.scatter(time_vector, recovered_count, label='Recuperados', color='green')

# Título y ejes
title = f"Modelo SIR discreto de red (N = {N}, b = {b}, 1/k = {1/k:.2f}, dt = {dt:.2f})"
plt.title(title)
plt.xlabel('Tiempo (días)')
plt.ylabel('Proporción de la población')

# Ajustar el eje x para marcas 
interval = 10  # Intervalo de 20 días para los ticks en el eje x
plt.xticks(ticks=np.arange(0, tmax + interval, interval))  # Ajustar para que se muestren los ticks cada 10 días
plt.xlim(0, tmax)  # Establecer el límite del eje x 
plt.legend()

# Agregar una cuadrícula
plt.grid(True)

plt.show()

#print('infected_neighbors_count: ', infected_neighbors_count)
#print('len(inf_n_c): ', len(infected_neighbors_count), 'tmax/dt: ', int(tmax/dt))

# Guardamos los vectores con el conteo de infectados, sanos y recuperados
base_name = f'red_b_{b}_k_{k:.2f}_N_{N}_dt_{dt:.2f}_neigh_{neighbors}'
np.savetxt(os.path.join(red_dir, f'susceptible_{base_name}.txt'), susceptible_count, fmt='%f', delimiter=',')
np.savetxt(os.path.join(red_dir, f'infected_{base_name}.txt'), infected_count, fmt='%f', delimiter=',')
np.savetxt(os.path.join(red_dir, f'recovered_{base_name}.txt'), recovered_count, fmt='%f', delimiter=',')
np.savetxt(os.path.join(red_dir, f'infected_neighbors_count_{base_name}.txt'), infected_neighbors_count, fmt='%f', delimiter=',')


## COMPARACIÓN CON EL MODELO DETERMINISTA (EDOS)
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
S, I, R = ret.T # noqa: E741

# Gráfica de resultados
plt.figure(figsize=(10, 6))
plt.plot(t, S, 'b', label='Susceptibles (Determinista)')
plt.plot(t, I, 'r', label='Infectados (Determinista)')
plt.plot(t, R, 'g', label='Recuperados (Determinista)')
title = f"Modelo SIR determinista vs discreto de red (# vecinos = {neighbors}, N = {N}, b = {b}, 1/k = {int(1/k)}, dt = {int(dt)})"
plt.title(title)
plt.xlabel('Tiempo (días)')
plt.ylabel('Proporción de la población')
plt.legend()  # Esto muestra la leyenda con las etiquetas definidas en 'label'

# Si queremos graficar solamente los puntos
plt.scatter(time_vector, susceptible_count, label='Susceptibles (Discreto - 4 vecinos)', color='blue', alpha=0.6)
plt.scatter(time_vector, infected_count, label='Infectados (Discreto - 4 vecinos)', color='red', alpha=0.6)
plt.scatter(time_vector, recovered_count, label='Recuperados(Discreto - 4 vecinos)', color='green', alpha=0.6)
interval = 10  # Intervalo de 20 días para los ticks en el eje x
plt.xticks(ticks=np.arange(0, tmax + interval, interval))  # Ajustar para que se muestren los ticks cada 10 días
plt.legend()  # Esto muestra la leyenda con las etiquetas definidas en 'label'
plt.xlim(0, tmax)  # Establecer el límite del eje x 
plt.show()