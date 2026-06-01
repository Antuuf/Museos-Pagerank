import numpy as np
from scipy.linalg import solve_triangular

def construye_adyacencia(D,m):
    # Función que construye la matriz de adyacencia del grafo de museos
    # D matriz de distancias, m cantidad de links por nodo
    # Retorna la matriz de adyacencia como un numpy.
    D = D.copy()
    l = [] # Lista para guardar las filas
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        l.append(fila<=fila[np.argsort(fila)[m]] ) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return(A)


def calculaLU(M):
    M = M.copy().astype(float)
    N = M.shape[0]

    U = M.copy()
    L = np.eye(N)

    for i in range(N):
        for j in range(i + 1, N):
            coeficiente = U[j, i] / U[i, i]
            U[j, :] -= coeficiente * U[i, :]
            L[j, i] = coeficiente

    return L, U

def calcula_matriz_C(A):
    # A: matriz de adyacencia (NxN)
    # Retorna: matriz de transiciones C = A^T @ K^-1 sin invertir explícitamente
    A_T = A.T
    suma_filas = np.sum(A, axis=1)  # Vector con la suma de cada fila de A. Representa a la diagonal K
    N = A.shape[0]

    C = np.zeros((N, N))

    for i in range(N):
        C[:, i] = A_T[:, i] / suma_filas[i]

    return C



def calcula_pagerank(A, alfa):
    N = A.shape[0]  # Cantidad de museos/nodos

    # 1. Calcular matriz de transiciones C
    C = calcula_matriz_C(A)

    # 2. Construir matriz M = (N / alfa) * (I - (1 - alfa) * C)
    I = np.eye(N)
    M = (N / alfa) * (I - (1 - alfa) * C)

    # 3. Construir vector b = vector de unos
    b = np.ones(N)

    # 4. Resolver Mp = b usando factorización LU
    L, U = calculaLU(M)

    # 5. Resolución por sustitución hacia adelante y atrás
    y = solve_triangular(L, b, lower=True)   # Ly = b
    p = solve_triangular(U, y, lower=False)               # Up = y

    return p

def calcula_matriz_C_continua(D):
    # Función para calcular la matriz de transiciones C en versión continua
    # D: Matriz de distancias entre museos

    D = D.copy()
    F = 1 / D
    np.fill_diagonal(F, 0)  # Evitamos dividir por 0 en la diagonal

    # K es la matriz diagonal con la suma de cada columna de F (según ecuación 4)
    suma_col = np.sum(F, axis=0)
    Kinv = np.diag(1 / suma_col)  # Inversa de la matriz K

    # Producto matricial: C = F @ Kinv (normaliza por columnas)
    C = F @ Kinv

    return C



def calcula_B(C, r):
    n = C.shape[0]
    B = np.eye(n)  # Empezamos con C^0 = I
    pot = np.eye(n)  # potencia acumulada de C

    for _ in range(1, r):
        pot = pot @ C
        B += pot

    return B


def norma_1_matriz(A):
    n = A.shape[1]
    max_suma = 0
    for j in range(n):
        suma_col = 0
        for i in range(A.shape[0]):
            suma_col += abs(A[i, j])
        if suma_col > max_suma:
            max_suma = suma_col
    return max_suma


def condicion_1_por_LU(B):
    norma_B = norma_1_matriz(B)

    # Factorización LU
    L, U = calculaLU(B)
    n = B.shape[0]
    norma_inv = 0

    # Resolver Bx = e_j para cada j
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1
        y = solve_triangular(L, e, lower=True)
        x = solve_triangular(U, y)

        norma_x = np.sum(np.abs(x))  # norma 1 de la columna de B^-1
        norma_inv = max(norma_inv, norma_x)

    return norma_B * norma_inv