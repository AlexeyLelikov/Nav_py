import numpy as np

def create_F(g : float,a : np.array, H : float, Lat : float, beta : float, T : float):
    R = 6356863 # м
    F = np.array([[1, 0, 0, -g * T, a[1] * T, 0, 0, 0],
                  [0, 1, g * T, 0, -a[0] * T, 0, 0, 0],
                  [0, -T / R, 1, 0, 0, T * np.cos(H), T * np.sin(H), 0],
                  [T / R, 0, 0, 1, 0, -T * np.sin(H), T * np.cos(H), 0],
                  [T * np.tan(Lat) / R, 0, 0, 0, 1, 0, 0, T],
                  [0, 0, 0, 0, 0, 1 - beta * T, 0, 0],
                  [0, 0, 0, 0, 0, 0, 1 - beta * T, 0],
                  [0, 0, 0, 0, 0, 0, 0, 1 - beta * T]])
    return F

def create_G(A : float, beta : float, T : float):
    filtr = A * np.sqrt(2 * beta) * T
    return np.diag([1,1,1,1,1,filtr,filtr,filtr])