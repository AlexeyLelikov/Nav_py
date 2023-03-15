import numpy as np
import pandas as pd
from scipy.io import loadmat
from Savage import Savage
import matplotlib.pyplot as plt
from LatLonAlt2XYZ import LatLonAlt2XYZ
from CalcGravity import CalcGravity
from SkewSymmMatr import SkewSymmMatr
from cross import cross
from Calc_F_C_N import Calc_F_C_N
from Calc_C_N_E import Calc_C_N_E
from twoPiBound import twoPiBound
from LatLonAlt2XYZ import LatLonAlt2XYZ
from DCM_bn import DCM_bn
from CalcDCMnue_align import CalcDCMnue_align
from EulerFromDCMnue import EulerFromDCMnue

mat = loadmat('nav.binB_03.mat')
data_df = pd.DataFrame(mat['UAV'])
data = data_df.to_numpy().T
data = np.delete(data, 0, 1)
data = data[405500:,:]
Sensors = np.zeros((6,2),dtype = np.float64)
n = data.shape[0]
acc = data[:,9:12]
gyro = data[:,12:15]

GPSQuality = data[:,23]
GPS = data[:,17:23]
GPS[:,0:2] = GPS[:,0:2] * np.pi / 180

StatusMS = data[:, 33]
MagnSens = data[:, 30:33] * 10
MagnSens[:,1:3] = -1 * MagnSens[:,1:3]

dt = 1 / 100

# ==================== Выставка ============================
Lat = GPS[60000,0]
Lon = GPS[60000,1]
Alt = GPS[60000,2]
Alignment = 60000
W_NUE = GPS[60000,3:6]
W_NUE_old = W_NUE
Lat_old = Lat
Alt_old = Alt
Roll = 0
Pitch = 0
Heading = 0
C_B_N = np.eye(3)
NavState = np.zeros(26)
MS_mean_count = 0
gyros_mean_count = 0
accels_mean = np.zeros(3)
gyros_mean = np.zeros(3)
MS_mean = np.zeros(3)
# ==================== Накопление данных =============================
for i in range(Alignment):
    accels_mean = accels_mean + acc[i,:]
    gyros_mean = gyros_mean + gyro[i,:]
    gyros_mean_count = gyros_mean_count + 1
    if (StatusMS[i] == 2):
        MS_mean = MS_mean + MagnSens[i]
        MS_mean_count = MS_mean_count + 1
# ======================== Расчет углов крена и тангажа ===============
accels_mean = (accels_mean/gyros_mean_count) / dt # Осреднение показаний акселерометров. Используется счетчик данных гироскопа
gyros_mean = (gyros_mean/gyros_mean_count) / dt # Осреднение показаний гироскопов + пересчет в рад/сек
MS_mean = (MS_mean / MS_mean_count)
if True: # если датчики - ВОГ, то проводится определение курса методом гирокомпасирования
    Rbn = CalcDCMnue_align([Lat,Lon,Alt],accels_mean,gyros_mean)
    (Heading, Pitch, Roll) = EulerFromDCMnue(Rbn)
    C_B_N = DCM_bn(Heading,Pitch,Roll) # Расчет DCM
    Heading = twoPiBound(Heading)  # Ограничение
# Вектор состояния навигационного алгоритма
NavState = np.array([Lat,Lon,Alt,          # Широта, долгота, геодезическая высота
                   W_NUE[0],W_NUE[1],W_NUE[2],               # Скорость на север, вверх, на восток
                   Roll,Pitch,Heading,   # Крен, тангаж, курс
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,           # Ошибки масштабных коэффициентов акселерометров, смещения нуля акселерометров, ошибки масштабных коэффициентов гироскопов
                   Lat_old,Alt_old,              # Предыдущие широта и высота
                   W_NUE_old[0],W_NUE_old[1],W_NUE_old[2]])
print("Параметры после выставки")
print("Курс ")
print(Heading * (180 / np.pi))
print("Тангаж ")
print(Pitch * (180 / np.pi))
print("Крен ")
print(Roll * (180 / np.pi))
print("Широта ")
print(Lat * (180 / np.pi))
print("Долгота ")
print(Lon * (180 / np.pi))
print("Высота ")
print(Alt)
print("W_N ")
print(W_NUE[0])
print("W_U ")
print(W_NUE[1])
print("W_E ")
print(W_NUE[2])
queue_Lat = np.zeros(n-(Alignment))
queue_Lon = np.zeros(n-(Alignment))
queue_pitch = np.zeros(n-(Alignment))
queue_roll = np.zeros(n-(Alignment))
queue_heading = np.zeros(n-(Alignment))
queue_Alt = np.zeros(n-(Alignment))
queue_W_N = np.zeros(n-(Alignment))
queue_W_U = np.zeros(n-(Alignment))
queue_W_E = np.zeros(n-(Alignment))
## ======================= Работа ============================
for i in range(Alignment,n):
    Sensors[0:3,1] = acc[i,:] # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
    Sensors[0:3,0] = acc[i-1,:]
    Sensors[3:6,0] = gyro[i-1,:] # выходной сигнал гироскопов (скорость) по оси (x,y,z) на предыдущем такте расчетов
    Sensors[3:6,1] = gyro[i,:] # выходной сигнал гироскопов (скорость) по оси (x,y,z) на текущем такте расчетов
    NavState[2] = 0
    NavState[4] = 0
    queue_Lat[i - Alignment] = NavState[0]
    queue_Lon[i - Alignment] = NavState[1]
    queue_roll[i - Alignment] = NavState[6]
    queue_pitch[i - Alignment] = NavState[7]
    queue_heading[i - Alignment] = NavState[8]
    queue_W_N[i - Alignment] = NavState[3]
    queue_W_U[i - Alignment] = NavState[4]
    queue_W_E[i - Alignment] = NavState[5]
    queue_Alt[i - Alignment] = NavState[2]
    (NavState, delta_v, delta_alpha, gt, C_B_N) = Savage(NavState, Sensors, dt, C_B_N)

n = queue_Lat.shape[0]
s = 60000
d = 230000
N = n-d
N_GPS = n-d-s
T_GPS = N_GPS / 100
T = N / 100 # c
t = np.linspace(0,T,N)
t_GPS = np.linspace(0,T_GPS,N_GPS)

queue_Lat = queue_Lat * (180 / np.pi)
queue_Lon = queue_Lon * (180 / np.pi)
queue_pitch = queue_pitch * (180 / np.pi)
queue_roll = queue_roll * (180 / np.pi)
queue_heading = queue_heading * (180 / np.pi)
queue_GPS_Lat = GPS[s:n-d,1]  * (180 / np.pi)
queue_GPS_Lon = GPS[s:n-d,2]  * (180 / np.pi)

plt.figure(1)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_Lat[:N])
plt.plot(t,queue_GPS_Lat[s:n-d])
plt.title("Широта")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(2)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_Lon[:N])
#plt.plot(t_GPS,queue_GPS_Lon)
plt.title("Долгота")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(3)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_pitch[:N])
plt.title("Тангаж")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(4)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_roll[:N])
plt.title("Крен")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(5)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_heading[:N])
plt.title("Курс")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(6)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_W_E[:N])
plt.title("Скорость N")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()

plt.figure(7)
plt.plot(t,queue_W_N[:N])
plt.title("Скорость E")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()

plt.figure(8)
fig = plt.gcf()
fig.set_size_inches(24, 10)
plt.plot(t,queue_W_U[:N])
plt.title("Скорость Up")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()

R = 6356.863 # км
queue_x = (queue_Lat * ((np.pi / 180) * R)) - (queue_Lat[0] * ((np.pi / 180) * R))
queue_y = (queue_Lon * (np.pi / 180) * (R * np.cos(Lat))) - (queue_Lon[0] * np.pi / 180 * (R * np.cos(Lat)))
plt.figure(9)
plt.plot(t,queue_x[:N])
plt.title("N")
plt.xlabel('Cекунды')
plt.ylabel('N,км')
plt.grid(True)
plt.show()

plt.figure(10)
plt.plot(t,queue_y[:N])
plt.title("E")
plt.xlabel('Cекунды')
plt.ylabel('E,км')
plt.grid(True)
plt.show()



