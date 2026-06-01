# Álgebra Lineal Computacional — Trabajo Práctico 1
## PageRank y Modelado de Flujo de Visitantes en la Noche de los Museos (CABA)

Este repositorio contiene la implementación y el análisis correspondientes al **Trabajo Práctico 1** de la materia **Álgebra Lineal Computacional (ALC)**. El proyecto se enfoca en el estudio de la red de museos de la Ciudad Autónoma de Buenos Aires (CABA) aplicando el algoritmo **PageRank** y modelos de flujo de visitantes basados en cadenas de Markov de tiempo discreto.

A partir de datos geográficos (coordenadas de latitud y longitud) y registros reales de visitas acumuladas, se analizan la centralidad de los museos y la dinámica de circulación del público durante el evento cultural "La Noche de los Museos".

---

##  Estructura del Proyecto

El repositorio está integrado por los siguientes archivos clave:

*   **[`TP1_template.ipynb`](file:///c:/Users/Antonella/Proyectos/ALC/TP1_template.ipynb)**: Notebook de Jupyter principal. Realiza la carga de datos geográficos (GeoJSON de museos y barrios de CABA), la visualización interactiva de mapas mediante `geopandas` y `matplotlib`, el análisis de grafos con `networkx`, y la ejecución de los experimentos numéricos y la resolución de sistemas lineales.
*   **[`template_funciones.py`](file:///c:/Users/Antonella/Proyectos/ALC/template_funciones.py)**: Módulo de Python con las implementaciones algorítmicas requeridas para el TP, optimizadas mediante operaciones vectorizadas en NumPy y resolución eficiente de sistemas triangulares usando SciPy.
*   **[`visitas.txt`](file:///c:/Users/Antonella/Proyectos/ALC/visitas.txt)**: Archivo con la cantidad de visitas registradas en cada uno de los museos evaluados.
*   **[`tp1.pdf`](file:///c:/Users/Antonella/Proyectos/ALC/tp1.pdf)**: Enunciado oficial detallado de la práctica con los fundamentos teóricos y las consignas de evaluación.

---

##  Modelado Matemático y Algoritmos

### 1. Construcción del Grafo y Matriz de Adyacencia ($A$)
Se genera un grafo dirigido donde cada nodo representa un museo. Para conectar los museos, se asume que un visitante en el museo $i$ solo elegirá desplazarse a uno de sus $m$ vecinos más cercanos (según distancia euclidiana). La matriz de adyacencia de tamaño $N \times N$ se define como:
$$A_{ij} = \begin{cases} 1 & \text{si } d(i, j) \le d(i, \text{m-ésimo vecino más cercano}) \text{ con } i \neq j \\ 0 & \text{en otro caso} \end{cases}$$
Se eliminan los autolazos (diagonal nula, $A_{ii} = 0$) para evitar transiciones triviales sobre el mismo museo.

### 2. Matriz de Transiciones ($C$)
*   **Modelo Discreto:** Las transiciones se distribuyen de manera uniforme entre las salidas disponibles:
    $$C = A^T K^{-1}$$
    Donde $K$ es la matriz diagonal que contiene los grados de salida (cantidad de enlaces salientes de cada museo).
*   **Modelo Continuo (Inversamente Proporcional a la Distancia):** La probabilidad de transitar a otro museo disminuye al aumentar la distancia física:
    $$F_{ij} = \begin{cases} \frac{1}{D_{ij}} & \text{si } i \neq j \\ 0 & \text{si } i = j \end{cases}$$
    La matriz de transiciones normalizada por columnas se obtiene como:
    $$C = F K^{-1} \quad \text{donde } K = \text{diag}(\sum_{i} F_{ij})$$

### 3. Sistema PageRank
El vector de centralidades (PageRank) $p$ se calcula al resolver el siguiente sistema de ecuaciones lineales:
$$M \cdot p = b$$
Donde:
*   $M = \frac{N}{\alpha} \left( I - (1 - \alpha) C \right)$
*   $b = \mathbf{1}$ (vector de unos de longitud $N$)
*   $\alpha \in (0, 1)$ es el factor de teletransportación (comportamiento del navegante aleatorio).
*   $I$ es la matriz identidad.

El sistema se resuelve aplicando **factorización LU** de la matriz de coeficientes $M$ ($M = L \cdot U$).

### 4. Flujo de Visitantes en la Noche de los Museos
Considerando que:
1.  $v$ es el vector que indica la cantidad de personas que ingresaron o comenzaron su recorrido en cada museo.
2.  Cada persona realiza exactamente $r$ visitas en total.
3.  $C$ es la matriz de transiciones.

El total de visitas acumuladas por cada museo al final de la noche se expresa como el vector $w$, relacionado mediante la matriz acumulada $B$:
$$w = v + C v + C^2 v + \dots + C^{r-1} v = B v \quad \text{con } B = I + C + C^2 + \dots + C^{r-1}$$
Para deducir la distribución inicial de personas ($v$) a partir de las visitas totales observadas ($w$), resolvemos:
$$B \cdot v = w$$
mediante factorización LU de $B$.

### 5. Análisis de Estabilidad y Cota del Error
La estabilidad numérica del sistema frente a perturbaciones en los datos medidos $w$ (por ejemplo, fallas en los sensores de conteo de ingresos) se analiza mediante el **número de condición** en norma 1 de la matriz $B$:
$$\text{cond}_1(B) = \|B\|_1 \cdot \|B^{-1}\|_1$$
Si se asume un error relativo conocido en $w$, la cota superior del error relativo propagado en el vector inicial recuperado $v$ cumple:
$$\frac{\|\Delta v\|_1}{\|v\|_1} \le \text{cond}_1(B) \cdot \frac{\|\Delta w\|_1}{\|w\|_1}$$

---

##  Biblioteca de Funciones (`template_funciones.py`)

El archivo contiene las siguientes rutinas desarrolladas desde cero:
*   `construye_adyacencia(D, m)`: Construye la matriz de adyacencia del grafo basándose en la matriz de distancias $D$ y el límite de enlaces $m$.
*   `calculaLU(M)`: Realiza la descomposición LU de una matriz cuadrada sin pivoteo.
*   `calcula_matriz_C(A)`: Computa la matriz de transiciones discreta a partir de la adyacencia.
*   `calcula_matriz_C_continua(D)`: Computa la matriz de transiciones continua usando la inversa de la distancia geográfica.
*   `calcula_pagerank(A, alfa)`: Construye la matriz del sistema PageRank y lo resuelve usando factorización LU y sustituciones triangulares.
*   `calcula_B(C, r)`: Calcula la sumatoria de potencias de transiciones $\sum_{k=0}^{r-1} C^k$.
*   `norma_1_matriz(A)`: Calcula de forma determinística la norma 1 de una matriz (máxima suma absoluta de columnas).
*   `condicion_1_por_LU(B)`: Calcula de forma analítica el número de condición en norma 1 de $B$, obteniendo su inversa mediante la resolución sistemática de sistemas de la forma $B x = e_j$ (para cada vector canónico $e_j$) mediante factorización LU.

---

##  Requisitos e Instrucciones de Uso

### Requisitos Previos
Es necesario contar con Python 3.8+ y las siguientes bibliotecas del ecosistema científico instaladas:

```bash
pip install numpy scipy pandas geopandas matplotlib seaborn networkx
```

*Nota: Para procesar correctamente el GeoJSON de barrios y museos de CABA dentro de `geopandas`, se requiere una conexión a internet activa durante la primera ejecución o descargar los archivos crudos en forma local.*

### Ejecución
Para visualizar el análisis y correr los experimentos numéricos de PageRank y flujo:

1. Clone el repositorio o colóquese en la carpeta del proyecto.
2. Inicie el entorno de Jupyter Notebook:
   ```bash
   jupyter notebook
   ```
3. Abra el archivo **`TP1_template.ipynb`** y corra las celdas secuencialmente.
