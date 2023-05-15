import numpy as np
from numpy import matmul as mn

def KF(X : np.array, # состояние (оцениваемые параметры)
       U : np.array, #  воздействия
       F : np.array, # матрица модели
       G : np.array, # матрица входного воздействия
       P : np.array, # матрица (коввариции) ошибки оценки
       Q : np.array, # матрица (коввариции) шума модели
       Z : np.array, # измерения
       H : np.array, # связь наблюдения и состояния
       R : np.array): # матрица (коввариции) шума измерения

    X_prior = mn(F,X) #+ mn(G,U)
    P_prior = mn(mn(F, P), F.T) + mn(mn(G, Q), G.T)

    S = mn(mn(H, P_prior), H.T) + R
    K = mn(mn(P_prior,H.T),np.linalg.inv(S))
    y = Z - mn(H, X_prior)

    X_posterior = X_prior + mn(K,y)
    P_posterior = mn((np.eye(len(X)) - mn(K,H)), P_prior)

    return X_posterior, P_posterior




