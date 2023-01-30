import numpy as np
from cross import cross
from SkewSymmMatr import SkewSymmMatr
from Calc_F_C_N import Calc_F_C_N
from twoPiBound import twoPiBound
# Статические переменные
Om_e = 7.292115e-5
E = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],dtype = np.float64)
# Используется ПЗ-90
a_e = 6378136.0
ecc = 0.0818191065283638
ecc2 = ecc * ecc
iter = 0
delta_v_array = np.zeros((3,4))
delta_alfa_array = np.zeros((3,4))

def SINS(NavState : np.array, Sensors : np.array, dt : float, C_i_b : np.array, C_i_n : np.array , gt : float) -> tuple:

    global Om_e
    global E
    # Используется ПЗ-90
    global a_e
    global ecc2
    global iter
    global delta_v_array
    global delta_alfa_array

    # переменные квази-скоростей и квази-угллов
    delta_v = Sensors[0:3] # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
    delta_alpha = Sensors[3:6] # выходной сигнал гироскопов (скорость) по оси (x,y,z) на текущем такте расчетов
    delta_v_array[0:3,iter] = delta_v
    delta_alfa_array[0:3,iter] = delta_alpha
    iter += 1

    if iter == 4:

        B = NavState[0]
        L = NavState[1]  # текущая  долгота ВС в выбранной СК
        h = NavState[2]  # текущая  геодезическая высота ВС относительно эллипсоида выбранной СК
        W = NavState[3:6] # текущий вектор земной скорости ВС (N,U,E)
        sinB = np.sin(B)
        cosB = np.cos(B)
        sinL = np.sin(L)
        cosL = np.cos(L)
        sinB2 = sinB * sinB
        R_N = a_e * (1 - ecc2) / (1 - ecc2 * sinB2) ** 1.5 + h
        R_E = a_e / np.sqrt(1 - ecc2 * sinB2) + h
        # Расчет вектора конечного поворота Эйлера (body)
        alfa_1_2 = np.array([delta_alfa_array[0,0],delta_alfa_array[1,0], delta_alfa_array[2,0]])+np.array([delta_alfa_array[0,1],delta_alfa_array[1,1], delta_alfa_array[2,1]])
        alfa_3_4 = np.array([delta_alfa_array[0,2],delta_alfa_array[1,2], delta_alfa_array[2,2]])+np.array([delta_alfa_array[0,3],delta_alfa_array[1,3], delta_alfa_array[2,3]])
        alfa = alfa_1_2 + alfa_3_4
        beta = 2 / 3 * cross(alfa_1_2,alfa_3_4) # коннинг
        alfa = alfa + beta # вектор порота body эйлера
        skew_alfa = SkewSymmMatr(alfa) # Коссосиметричная матрица
        mod_alpha = np.sqrt(alfa[0] * alfa[0] + alfa[1] * alfa[1] + alfa[2] * alfa[2])
        mod_alpha2 = mod_alpha * mod_alpha
        # Расчет вектора конечного поворота (скорости) Эйлера (Navigation)
        omega_IE_N = np.array([Om_e * cosB, Om_e * sinB, 0.0 ])
        F_C_N = Calc_F_C_N(B, h, a_e, ecc2)
        omega_sum = omega_IE_N + F_C_N @ W
        skew_omega = SkewSymmMatr(omega_sum)
        # нахождение МНК из body t-1 в body t при помощи решения уравенения Пуасона
        C_b_b = E + (1 - mod_alpha2 / 6) * skew_alfa + (0.5 - mod_alpha2/24) * skew_alfa @ skew_alfa
        C_i_b = C_b_b @ C_i_b  # Получили МНК из body в инерциальную в данный такт
        # нахождение МНК из навигационной СК t-1 в навигационную СК t при помощи решения уравенения Пуасона
        C_n_n = E - dt * skew_omega + (dt**2) / 2 * skew_omega @ skew_omega
        C_i_n = C_n_n @ C_i_n
        C_b_n = C_i_n @ C_i_b.T
        #нормализация и ортогонализация мнк
        #C_i_b = (E - ((1.0 / 2.0) * ( C_i_b @  C_i_b.T - E))) @  C_i_b
        #C_i_n = (E - ((1.0 / 2.0) * ( C_i_n @  C_i_n.T - E))) @  C_i_n
        #C_b_n = (E - ((1.0 / 2.0) * ( C_b_n @  C_b_n.T - E))) @  C_b_n
        #нахождение приращения линейных скоростей методом Рунге-Кутта 4-го порядка
        delta_W = np.array([0,0,0])
        for i in range(4):
            skew_delta_alfa = SkewSymmMatr(delta_alfa_array[0:3,i])
            w_b = delta_v_array[0:3,i]
            k_1 = w_b - skew_delta_alfa @ delta_W
            k_2 = w_b - skew_delta_alfa @ (delta_W + (dt / 2) * k_1)
            k_3 = w_b - skew_delta_alfa @ (delta_W + (dt / 2) * k_2)
            k_4 = w_b - skew_delta_alfa @ (delta_W + dt * k_3)
            delta_W = delta_W + (k_1+2*k_2+2*k_3+k_4)/6

        delta_W_n = C_b_n @ delta_W

        omega = omega_IE_N + F_C_N @ W
        W_N = W[0] + delta_W_n[0] + dt * (-(omega_IE_N[2] + omega[2]) * W[2] + omega[2] * W[1])
        W_U = W[1] + delta_W_n[1] + dt * ((omega_IE_N[0] + omega[0]) * W[2] + omega[2] * W[0] - gt)
        W_E = W[2] + delta_W_n[2] + dt * ((omega_IE_N[2] + omega[2]) * W[0] - (omega_IE_N[0] + omega[0])*W[1])

        Lat = B + dt * W_N / R_N # Широта
        Lon = L + dt * W_E / R_E # долгота
        Alt = 0
        roll = np.arctan2(-C_b_n[1,2], C_b_n[1, 1])
        pitch = np.arctan(C_b_n[1, 0] / np.sqrt(C_b_n[1, 1] * C_b_n[1, 1] + C_b_n[1, 2] * C_b_n[1, 2]))
        heading = np.arctan2(C_b_n[2, 0], C_b_n[0, 0])
        heading = twoPiBound(heading) # Ограничение

        NavState = np.array([Lat,Lon,Alt,          # Широта, долгота, геодезическая высота
                   W_N,W_U,W_E,            # Скорость на север, вверх, на восток
                    roll,pitch,heading])           # Крен , тангаж , курс

        delta_v_array = np.zeros((3,4))
        delta_alfa_array = np.zeros((3,4))
        iter = 0

    return NavState , C_i_b , C_i_n