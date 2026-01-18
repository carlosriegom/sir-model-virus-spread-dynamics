# Modelo SIR: Dinámica de propagación de un virus

**Autores:** Pablo Rodríguez Soria y José Carlos Riego

Este proyecto estudia la dinámica de propagación de un virus utilizando el modelo epidemiológico **SIR** (Susceptibles-Infectados-Recuperados) desde la perspectiva de la **Mecánica Estadística**. El trabajo compara soluciones deterministas (ecuaciones diferenciales) con modelos estocásticos (Monte Carlo) en campo medio y en redes discretas, analizando finalmente el impacto de diferentes estrategias de confinamiento.

## Tabla de Contenidos

1. [Fundamentos Teóricos: Modelo SIR Determinista](#1-fundamentos-teóricos-modelo-sir-determinista)
2. [Modelo de Campo Medio](#2-modelo-de-campo-medio)
3. [Modelo Discreto de Red](#3-modelo-discreto-de-red)
4. [Escenarios de Confinamiento](#4-escenarios-de-confinamiento)
5. [Instalación y Uso](#5-instalación-y-uso)

---

## 1. Fundamentos Teóricos: Modelo SIR Determinista

El modelo base es un sistema determinista descrito por Ecuaciones Diferenciales Ordinarias (EDOs) que rigen la dinámica no hamiltoniana de un sistema cerrado. La población se divide en tres compartimentos: **S** (Susceptibles), **I** (Infectados) y **R** (Recuperados).

### Ecuaciones del sistema

La evolución temporal de las proporciones de población está dada por:

$$ \frac{ds(t)}{dt} = -b \cdot i(t) \cdot s(t) $$

$$ \frac{di(t)}{dt} = b \cdot i(t) \cdot s(t) - k \cdot i(t) $$

$$ \frac{dr(t)}{dt} = k \cdot i(t) $$

### Donde:

- $b \in (0, 1]$: **Tasa de transmisión** del virus ($días^{-1}$).
- $k$: Tasa de recuperación, definida como el inverso del tiempo promedio de recuperación $t_p = 1/k$.

**Resultados clave:** Las simulaciones muestran cómo variar $b$ y $t_p$ afecta la velocidad de la infección y el pico de la curva de infectados.

---

## 2. Modelo de Campo Medio

Se implementa una aproximación estocástica utilizando el **Método de Monte Carlo**.

- **Inicialización:** Matriz `network` donde cada individuo tiene un estado (1: Susceptible, 2: Infectado, 3: Recuperado).
- **Probabilidades de transición:**
  - **Infección:** Depende de la proporción global de infectados $i(t)$. La probabilidad es $P_i(t) = b \cdot \Delta t \cdot i(t)$.
  - **Recuperación:** Se modela con un tiempo característico $\tau$. La probabilidad es $P_r(t) = 1 - e^{-t/\tau}$.

**Conclusión:** Para poblaciones grandes ($N \gg 1$), el comportamiento del modelo de campo medio converge y se aproxima fielmente a la solución analítica del modelo determinista.

---

## 3. Modelo Discreto de Red

Este modelo introduce la **topología espacial**. Los individuos están fijos en una red y solo interactúan con sus vecinos directos, simulando el contacto local.

### Dinámica Local

La probabilidad de contagio ya no depende del promedio global, sino de los vecinos locales infectados ($v$):

$$ P_c = 1 - (1 - p)^v $$

Donde $p$ es la probabilidad de contagio en un solo contacto.

### Calibración de Parámetros ($b$ vs $p$)

Para comparar el modelo de red con el determinista, se estableció una relación empírica entre la tasa de transmisión global $b$ y la probabilidad local $p$. Mediante ajuste logarítmico de los datos simulados para 4 vecinos, se obtuvo:

$$ p(b) = 0.21 \cdot \log(b) + 0.71 $$

Esto permite traducir parámetros macroscópicos a comportamientos locales.

### Topologías

Se analizaron redes de **4 vecinos** y **8 vecinos**. Aumentar el número de vecinos (contactos) acelera drásticamente la curva de contagios.

---

## 4. Escenarios de Confinamiento

Se simularon estrategias de distanciamiento social modificando la conectividad de la red a partir de un día específico ($t = 20$).

| Tipo de Confinamiento | Descripción               | Efecto en la Red                                                           |
| :-------------------- | :------------------------ | :------------------------------------------------------------------------- |
| **Permisivo**         | Restricción por columnas. | Los individuos interactúan con **2 vecinos** (arriba/abajo).               |
| **Estricto**          | Restricción por parejas.  | Los individuos interactúan con **1 vecino** (pareja dentro de la columna). |

### Resultados:

- Ambos confinamientos logran "aplanar la curva".
- El confinamiento **estricto** reduce la proporción de infectados de manera mucho más eficiente y rápida que el permisivo, deteniendo la propagación casi por completo poco después de su implementación.

---

## 5. Instalación y Uso

Para ejecutar los diferentes scripts del proyecto, es necesario tener `Conda`o `Miniconda` instalado, y crear un entorno virtual:

```bash
conda env create -f environment.yml
```
