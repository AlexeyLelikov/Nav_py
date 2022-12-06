import numpy as np
import math
from SkewSymmMatr import SkewSymmMatr
from cross import cross
from LatLonAlt2XYZ import LatLonAlt2XYZ
from CalcGravity import CalcGravity
from Calc_F_C_N import Calc_F_C_N
from Calc_C_N_E import Calc_C_N_E
from twoPiBound import twoPiBound

def Savage(NavState : np.array, Sensors : np.array , dt : float  , C_B_N_old : np.array ) -> tuple:

    Om_e = 7.292115e-5
    EYE3x3 = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = np.float64)

    # Используется ПЗ-90
    a_e = 6378136
    ecc = 0.0818191065283638
    ecc2 = ecc * ecc

    B = NavState[0]  # текущая  широта ВС в выбранной СК
    L = NavState[1]  # текущая  долгота ВС в выбранной СК
    h = NavState[2]  # текущая  геодезическая высота ВС относительно эллипсоида выбранной СК

    B_old = NavState[21]  # предыдущее значение широты ВС
    h_old = NavState[22]  # предыдущее значение геодезической высоты ВС

    W = NavState[3:6] # текущий вектор земной скорости ВС (N,U,E)
    W_old = NavState[23:26] # предыдущий вектор земной скорости ВС (N,U,E)

    # переменные квази-скоростей и квази-угллов
    delta_v_old = Sensors[0:3,0] #  выходной сигнал акселерометра по осям чувствительности (x,y,z) на предыдущем такте расчетов
    delta_v = Sensors[0:3,1] # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
    delta_alpha_old = Sensors[3:6,0] # выходной сигнал гироскопов (скорость) по оси (x,y,z) на предыдущем такте расчетов
    delta_alpha = Sensors[3:6,1] #  выходной сигнал гироскопов (скорость) по оси (x,y,z) на текущем такте расчетов

    # Скалинг
    mod_alpha = np.sqrt(delta_alpha[0] * delta_alpha[0] + delta_alpha[1] * delta_alpha[1] + delta_alpha[2] * delta_alpha[2])
    mod_alpha2 = mod_alpha * mod_alpha

    beta = 1 / 12 * cross(delta_alpha_old, delta_alpha)

    Sa = (dt / 12) * (5 * delta_alpha + delta_alpha_old)
    Sv = (dt / 12) * (5 * delta_v + delta_v_old)

    temp1 = delta_alpha - delta_alpha_old
    temp2 = delta_v - delta_v_old

    dR_scrl_A =  - 1 / 24 * ((cross(temp1,Sv) + cross(temp2,Sa)))
    dR_scrl_B = dt / 144 * cross(temp2,delta_alpha) - dt / 144 * cross(temp1,delta_v) + dt / 240 * cross(temp1,temp2)
    delta_R_scrl = dR_scrl_A + dR_scrl_B

    B_extrap = 1.5 * B - 0.5 * B_old
    h_extrap = 1.5 * h - 0.5 * h_old
    W_extrap = 1.5 * W - 0.5 * W_old

    sinB = np.sin(B_extrap)
    cosB = np.cos(B_extrap)
    sinL = np.sin(L)
    cosL = np.cos(L)

    sinB2 = sinB * sinB

    gt = 9.7803267715 * (1 + 0.0052790414 * sinB2 + 0.0000232718 * sinB2 * sinB2) + (-0.0000030876910891 + 0.0000000043977311 * sinB2) * h_extrap + 0.0000000000007211 * h_extrap * h_extrap

    arr_X_Y_Z = LatLonAlt2XYZ(B_extrap, L, h_extrap)
    X = arr_X_Y_Z[0]
    Y = arr_X_Y_Z[1]
    Z = arr_X_Y_Z[2]

    gravity = CalcGravity(X, Y, Z)

    g_P_N_extrap = np.zeros((3))

    g_P_N_extrap[0] = -cosL * sinB * gravity[0] - sinL * sinB * gravity[1] + cosB * gravity[2]
    g_P_N_extrap[1]= cosL * cosB * gravity[0] + sinL * cosB * gravity[1] + sinB * gravity[2]

    gt = abs(g_P_N_extrap[1])

    omega_IE_N_extrap = np.array([Om_e * cosB, Om_e * sinB, 0 ])

    F_C_N_extrap = Calc_F_C_N(B_extrap, h_extrap, a_e, ecc2)

    if mod_alpha > 1e-16:
        temp1 = (1 - np.cos(mod_alpha)) / mod_alpha2
        temp2 = 1 / mod_alpha2 * (1 - np.sin(mod_alpha) / mod_alpha)
    else:
        temp1 = 1 / 2
        temp2 = 1 / 6

    delta_v_rot_m = temp1 * cross(delta_alpha, delta_v) + temp2 * cross(delta_alpha, cross(delta_alpha, delta_v))

    delta_v_scul = 1 / 12 * cross(delta_alpha_old, delta_v) + 1 / 12 * cross(delta_v_old, delta_alpha)

    delta_v_SFm_BmMinus1 = delta_v + delta_v_rot_m + delta_v_scul

    delta_v_SFm_LnMinus1 = C_B_N_old @ delta_v_SFm_BmMinus1

    ksi_nMinus1_m = omega_IE_N_extrap * dt + (F_C_N_extrap @ W_extrap) * dt

    C_LnMinus1_Lm = EYE3x3 - SkewSymmMatr(ksi_nMinus1_m)

    delta_v_SFm_L = C_LnMinus1_Lm @ delta_v_SFm_LnMinus1

    delta_v_G_CORm_N = g_P_N_extrap * dt - cross((2 * omega_IE_N_extrap + F_C_N_extrap @ W_extrap), W_extrap) * dt

    W_new = W + delta_v_SFm_L + delta_v_G_CORm_N

    W_new[1] =  W_new[1] * 0.9999999

    temp3 = cross(Sa, delta_v) + cross(delta_alpha, Sv)

    delta_R_Rotm = 1 / 6 * temp3

    delta_R_SFm_B = Sv + delta_R_Rotm + delta_R_scrl

    delta_R_SFm_L = - 1 / 6 * cross(ksi_nMinus1_m, delta_v_SFm_LnMinus1) * dt + C_B_N_old @ delta_R_SFm_B

    delta_R_m_N = (W + 1 / 2 * delta_v_G_CORm_N) * dt + delta_R_SFm_L

    delta_h = delta_R_m_N[1]

    h_new = h + delta_h

    ksi_n = F_C_N_extrap @ delta_R_m_N

    mod_ksi_n = np.sqrt(ksi_n[0] * ksi_n[0] + ksi_n[1] * ksi_n[1] + ksi_n[2] * ksi_n[2])

    mod_ksi_n_2 = mod_ksi_n * mod_ksi_n

    if mod_ksi_n > 1e-16:
        temp1 = np.sin(mod_ksi_n) / mod_ksi_n
        temp2 = (1 - np.cos(mod_ksi_n)) / mod_ksi_n_2
    else:
        temp1 = 1
        temp2 = 1 / 2

    C_Nn_NnMinus1 = EYE3x3 + temp1 * SkewSymmMatr(ksi_n) + temp2 * SkewSymmMatr(ksi_n) @ SkewSymmMatr(ksi_n)

    C_N_E_old = Calc_C_N_E(B, L)

    C_N_E_new = C_N_E_old @ C_Nn_NnMinus1

    C_N_E_new = (EYE3x3 - 1 / 2 * (C_N_E_new @ C_N_E_new.T - EYE3x3)) @ C_N_E_new

    B_new = np.atan2(C_N_E_new[2, 1], C_N_E_new[2, 0])

    L_new = np.atan2(- C_N_E_new[0, 2], C_N_E_new[1, 2])

    phi_m = delta_alpha + beta

    mod_phi_m = np.sqrt(phi_m[0] * phi_m[0]+ phi_m[1] * phi_m[1] + phi_m[2] * phi_m[2])

    mod_phi_m_2 = mod_phi_m * mod_phi_m

    if mod_phi_m > 1e-16:
        temp1 = np.sin(mod_phi_m) / mod_phi_m
        temp2 = (1 - np.cos(mod_phi_m)) / mod_phi_m_2
    else:
        temp1 = 1
        temp2 = 1 / 2

    C_Bm_BmMinus1 = EYE3x3 + temp1 * SkewSymmMatr(phi_m) + temp2 * SkewSymmMatr(phi_m) @ SkewSymmMatr(phi_m)

    C_Bm_LnMinus1 = C_B_N_old @ C_Bm_BmMinus1

    B_mean = 1 / 2 * (B_new + B)
    h_mean = 1 / 2 * (h_new + h)

    omega_IE_N_mean = np.array([Om_e * np.cos(B_mean), Om_e * np.sin(B_mean), 0])

    F_C_N_mean = Calc_F_C_N(B_mean, h_mean, a_e, ecc2)

    zetta_n = omega_IE_N_mean * dt + F_C_N_mean @ delta_R_m_N

    mod_zetta_n = np.sqrt(zetta_n[0] * zetta_n[0] + zetta_n[1] * zetta_n[1] + zetta_n[2] * zetta_n[2])

    mod_zetta_n_2 = mod_zetta_n * mod_zetta_n

    if mod_zetta_n > 1e-16:
        temp1 = np.sin(mod_zetta_n) / mod_zetta_n
        temp2 = (1 - np.cos(mod_zetta_n)) / mod_zetta_n_2
    else:
        temp1 = 1
        temp2 = 1 / 2

    C_LnMinus1_Ln = EYE3x3 - temp1 * SkewSymmMatr(zetta_n) + temp2 * SkewSymmMatr(zetta_n) * SkewSymmMatr(zetta_n)

    C_Bm_Ln = C_LnMinus1_Ln @ C_Bm_LnMinus1

    C_B_L_new = C_Bm_Ln

    roll = np.atan2(-C_B_L_new[1,2], C_B_L_new[1, 1])
    pith = np.atan(C_B_L_new[1, 0] / np.sqrt(C_B_L_new[1, 1] * C_B_L_new[1, 1] + C_B_L_new[1, 2] * C_B_L_new[1, 2]))
    heading = np.atan2(C_B_L_new[2, 0], C_B_L_new[0, 0])
    heading = twoPiBound(heading) # Ограничение

    NavState[0] = B_new
    NavState[21] = B
    NavState[1] = L_new
    NavState[2] = h_new
    NavState[22] = h
    #NavState[3:6] = W_new
    NavState[3] = W_new[0]
    NavState[4] = W_new[1]
    NavState[5] = W_new[2]
    #NavState[23:26] = W
    NavState[23] = W[0]
    NavState[24] = W[0]
    NavState[25] = W[0]
    NavState[6] = roll
    NavState[7] = pith
    NavState[8] = heading

    return NavState , delta_v, delta_alpha, gt , C_B_N

Lat = 55 * (np.pi / 180)
Lon = 38 * (np.pi / 180)

Lat_old = Lat
Alt = 0
Alt_old = 0

# Вычисление силы тяжести

arr_X_Y_Z = LatLonAlt2XYZ(Lat, Lon, Alt)
X = arr_X_Y_Z[0]
Y = arr_X_Y_Z[1]
Z = arr_X_Y_Z[2]
gravity = CalcGravity(X, Y, Z)
sinB = np.sin(Lat)
cosB = np.cos(Lat)
sinL = np.sin(Lon)
cosL = np.cos(Lon)

g_P_N_extrap = np.zeros(3)
g_P_N_extrap[0] = -cosL * sinB * gravity[0] - sinL * sinB * gravity[1] + cosB * gravity[2]
g_P_N_extrap[2] = cosL * cosB * gravity[0] + sinL * cosB * gravity[1]+ sinB * gravity[2]
gt = abs(g_P_N_extrap[1])
U = 7.292115e-5
acc = np.array([0.0, gt / 100.0, 0.0])
gyro = np.array([(U * np.cos(Lat)) / 100.0, (U * np.sin(Lat)) / 100.0, 0.0])
Sensors = np.zeros((6,2))
W_NUE = np.array([0,0,0])
W_NUE_old = np.array([0,0,0])
dt = 1.0 / 100.0
Roll = 0.0
Pitch = 0.0
Heading = 0
T = 360
N = 100 * T
t = np.linspace(0,T,N)
C_B_N = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = np.float64)

NavState = np.array([Lat,Lon,Alt,          # Широта, долгота, геодезическая высота
                   W_NUE[0],W_NUE[1],W_NUE[2],               # Скорость на север, вверх, на восток
                   Roll,Pitch,Heading,   # Крен, тангаж, курс
                   0,0,0,0,0,0,0,0,0,0,0,0,           # Ошибки масштабных коэффициентов акселерометров, смещения нуля акселерометров, ошибки масштабных коэффициентов гироскопов
                   Lat_old,Alt_old,              # Предыдущие широта и высота
                   W_NUE_old[0],W_NUE_old[1],W_NUE_old[2]])

queue_Lat = np.zeros(N)
queue_Lon = np.zeros(N)
queue_pitch = np.zeros(N)
queue_roll = np.zeros(N)
queue_yaw = np.zeros(N)
queue_Alt = np.zeros(N)
queue_W_N = np.zeros(N)
queue_W_U = np.zeros(N)
queue_W_E = np.zeros(N)

Sensors[0:3,1] = acc # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
Sensors[0:3,0] = acc
Sensors[3:6,0] = gyro
Sensors[3:6,1] = gyro

for i in range(0,N):
    NavState[2] = 0.0
    NavState[4] = 0.0
    queue_Lat[i] = NavState[0]
    queue_Lon[i] = NavState[1]
    queue_roll[i] = NavState[6]
    queue_pitch[i] = NavState[7]
    queue_yaw[i] = NavState[8]
    queue_W_N[i] = NavState[3]
    queue_W_U[i] = NavState[4]
    queue_W_E[i] = NavState[5]
    queue_Alt[i] = NavState[2]
    (NavState, delta_v, delta_alpha, gt, C_B_N) = Savage(NavState, Sensors, dt, C_B_N)




