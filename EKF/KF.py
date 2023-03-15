import numpy as np
from numpy import matmul as mn
from scipy import integrate
import matplotlib.pyplot as plt
import math
from math import cos , sin
from copy import copy

def KF(X = np.array([[]]), # состояние (оцениваемые параметры)
    U = np.array([[]]), #  воздействия
    F = np.array([[]]), # матрица модели
    G = np.array([[]]), # матрица входного воздействия
    P = np.array([[]]), #
    Q = np.array([[]]), # матрица (коввариции) шума модели
    Gamma = np.array([[]]), #
    Z = np.array([[]]), # наблюдение
    H = np.array([[]]), # связь наблюдения и состояния
    R = np.array([[]])):

    X_prior = mn(F,X) + mn(G,U)
    P_prior = mn(mn(F, P), F.T) + mn(mn(Gamma, Q), Gamma.T)

    y = Z - mn(H, X_prior)

    S = mn(mn(H, X_prior), H.T) + R
    K = mn(mn(P_prior,H.T),np.linalg.inv(S))

    X_posterior = X_prior + mn(K,y)
    P_posterior = mn((np.eye(len(X)) - mn(K,H)), P_prior)

    return  X_posterior , P_posterior

t = np.linspace(0,8,100)
delta_t = t[1] - t[0]

F = np.array([[1, 0, delta_t,       0],
              [0, 1,       0, delta_t],
              [0, 0,       1,       0],
              [0, 0,       0,       0]])

G = np.array([[delta_t**2,          0],
              [         0, delta_t**2],
              [delta_t**2,          0],
              [         0, delta_t**2]])

Q = np.eye(2)

def U_x(t):
    return 2.5 * sin(t)

def U_y(t):
    return -2*sin(t)

def se_solve(Y,t):
    x,y,v_x,x_y = Y
    return []




