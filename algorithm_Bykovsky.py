import matplotlib.pyplot as pyplot
import numpy as np
import pickle #для сохранения данных в bin файл 

deg_to_rad = np.pi/180
rad_to_deg = 1/deg_to_rad

# Загрузка из файла

with open('modelled_sensetive_elements.pkl', 'rb') as file:
     vars = pickle.load(file)

#задаю углы ориентации
#     #для работы
psi = vars['psi']#курс в рад
gamma = vars['gamma'] #крен
theta = vars['theta']#тангаж

#     # для экспериментов
# psi = 0 *deg_to_rad#курс
# gamma = 0*deg_to_rad #крен
# theta = 0*deg_to_rad #тангаж
#для красивого вывода графиков
psi_0 = psi
gamma_0= gamma
theta_0 = theta


# задаю начальные координаты (Москва, Оснанкинская башня)
phi_0 = vars['lattitude'] # в радианах
lmbda_0 = vars['longitude'] # в радианах

g = 9.81
U = 15 * deg_to_rad / 3600 #радианы в секунды
R_e = 6400*10**3 #м экваториальный радиус
a =	6378245 #большая полуось
b = 6356856
e = np.sqrt(1-b**2/a**2) #эксцентрисетет
R_lambda=R_e/((1-e**2*np.sin(phi_0)**2)**(1/2)) #радиус запад-восток
R_phi=R_e*(1-e**2)/((1-e**2*np.sin(phi_0)**2)**(3/2))#радиус север-юг

E = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) #единичная матрица

#компенсацию смещения нуей и ошибок мы пропускаем

#для гироскопов у нас выдаются углы приращения за один такт
freq = vars['freq'] #частота

h_1 = 1/freq #время такта сверхбыстрого цикла
h_4=4*h_1
h_6 = 6*h_1

delta_alpha_bx = np.array(vars['omega_x']) * h_1
delta_alpha_by = np.array(vars['omega_y'])*h_1
delta_alpha_bz = np.array(vars['omega_z'])*h_1

w_bx=np.array(vars['a_x'])*h_1
w_by=np.array(vars['a_y'])*h_1
w_bz=np.array(vars['a_z'])*h_1

delta_W = np.array([[0], [0], [0]])

V_ox = 0
V_oy = 0
V_oz =0 
Hight = 0 #высота (инерциальная)
Hight_ist = 135#истинная высота (баро)
c1, c2 = 4, 5#для вертикального канала
omega_ox, omega_oy, omega_oz = 0,0,0
C_i_o = E

C_b_o = np.array([[np.cos(gamma)*np.cos(psi)+np.sin(gamma)*np.sin(theta)*np.sin(psi), np.cos(theta)*np.sin(psi), np.sin(gamma)*np.cos(psi)-np.cos(gamma)*np.sin(psi)*np.sin(theta)],
[-np.cos(gamma)*np.sin(psi)+np.sin(gamma)*np.sin(theta)*np.cos(psi),np.cos(theta)*np.cos(psi), -(np.sin(gamma)*np.sin(psi)+np.cos(gamma)*np.sin(theta)*np.cos(psi))],
[-np.sin(gamma)*np.cos(theta), np.sin(theta), np.cos(gamma)*np.cos(theta)]])

C_o_b=C_b_o.T #транспонированная матрица

C_i_b = C_b_o.T

phi = phi_0
lmbda=lmbda_0

#массивы скоростей, координат и углов

V_x, V_y = [], []
phi_nav, lmbda_nav =[], []
psi, gamma, theta = [], [], []

sum_a_bx,sum_a_by, sum_a_bz = 0,0,0
flag=True
num = 0
for i in range(len(vars['omega_x'])-4):
# for i in range(5*60*100):
    #вычисление вектора Эйлера
    V_oz = 0
    Hight = 0
    if num >=4:
        flag = False
        num = 0
        sum_a_bx,sum_a_by, sum_a_bz = 0,0,0
    if flag:
        if num<4:
            sum_a_bx+=delta_alpha_bx[i]
            sum_a_by+=delta_alpha_by[i]
            sum_a_bz+=delta_alpha_bz[i]
            num+=1
            # continue
        # else:
        #     flag =False
    else:
    # sum_a_bx = np.sum(delta_alpha_bx[i:i+4])
    # sum_a_by = np.sum(delta_alpha_by[i:i+4])
    # sum_a_bz = np.sum(delta_alpha_bz[i:i+4])
    # A = 2/3*np.array([[(delta_alpha_by[i] + delta_alpha_by[i+1])*(delta_alpha_bz[i+3]+delta_alpha_bz[i+4]) - (delta_alpha_by[i+2] + delta_alpha_by[i+3])*(delta_alpha_bz[i]+delta_alpha_bz[i+1])],    [(delta_alpha_bx[i+2] + delta_alpha_bx[i+3])*(delta_alpha_bz[i]+delta_alpha_bz[i+1]) - (delta_alpha_bx[i]+delta_alpha_bx[i+1])*(delta_alpha_bz[i+2]+delta_alpha_bz[i+3])],    [(delta_alpha_bx[i] + delta_alpha_bx[i+1])*(delta_alpha_bz[i+2] + delta_alpha_bz[i+3]) - (delta_alpha_bx[i+2]+delta_alpha_bx[i+3])*(delta_alpha_bz[i] + delta_alpha_bz[i+1])]])
        vec1 = np.array([delta_alpha_bx[i],delta_alpha_by[i], delta_alpha_bz[i]])+np.array([delta_alpha_bx[i+1],delta_alpha_by[i+1], delta_alpha_bz[i+1]])
        vec2 = np.array([delta_alpha_bx[i+2],delta_alpha_by[i+2], delta_alpha_bz[i+2]])+np.array([delta_alpha_bx[i+3],delta_alpha_by[i+3], delta_alpha_bz[i+3]])
        A = 2/3*np.cross(vec1,vec2).reshape(3,1)
        Theta = np.array([[sum_a_bx],[sum_a_by],[sum_a_bz]])+ A #матрица столбец

        Theta_sym = np.array([[0, -Theta[2,0], Theta[1,0]], [Theta[2,0], 0, -Theta[0,0]], [-Theta[1,0], Theta[0,0], 0]])
        Theta_amp = np.sqrt(Theta[0,0]**2 + Theta[1,0]**2 + Theta[2,0]**2)
        #матрица
    
        C_b_b = E+ (1- Theta_amp**2/6)*Theta_sym + (0.5 - Theta_amp**2/24)*np.dot(Theta_sym, Theta_sym)

        #потом продолжим
        omega_sym = np.array([[0, -omega_oz, omega_oy], [omega_oz, 0, -omega_ox], [-omega_oy, omega_ox, 0]])

        delta_C_o_o = E - h_4*omega_sym + h_4**2/2*np.dot(omega_sym,omega_sym)

        C_i_o = np.dot(delta_C_o_o, C_i_o)
        C_b_i = np.dot(C_b_b, C_i_b)
        C_b_o = np.dot(C_i_o, C_b_i.T)

        #нормализация и ортогонализация матрицы

        #сначала нормализация для строк
        for i in range(3):
            delta_ii = 1 - np.dot(C_b_o[i,:], C_b_o[i,:].T)
            if abs(delta_ii)>10**(-5):
                # C_b_o[i,:] = C_b_o[i,:] - 0.5*delta_ij*C_b_o[i,:]
                C_b_o[i,:] = C_b_o[i,:]/np.sqrt(np.dot(C_b_o[i,:], C_b_o[i,:].T))
        #потом нормализация для столбцов
        for i in range(3):
            delta_ii = 1 - np.dot(C_b_o[:,i], C_b_o[:,i].T)
            if abs(delta_ii)>10**(-5):
                # C_b_o[:,i] = C_b_o[:,i] - 0.5*delta_ij*C_b_o[:,i]
                C_b_o[:,i] = C_b_o[:,i]/np.sqrt(np.dot(C_b_o[:,i], C_b_o[:,i].T))

        #ортогонализация
    
        #столбцы
        if (i%10 in [0,1,2,3,4,5]):
            for i in range(3):
                for j in range(3):
                    if i!=j:
                        delta_ij = np.dot(C_b_o[:,i], C_b_o[:,j])
                        if abs(delta_ij)>1e-14:
                            C_b_o[:,j] = C_b_o[:,j] - 0.5*delta_ij*C_b_o[:,i]
                            C_b_o[:,i] = C_b_o[:,i] - 0.5*delta_ij*C_b_o[:,j]
                            # print(np.dot(C_b_o[:,i], C_b_o[:,j]))
        else:
            #строки
            for i in range(3):
                for j in range(3):
                    if i!=j:
                        delta_ij = np.dot(C_b_o[i,:], C_b_o[j,:])
                        if abs(delta_ij)>1e-14:
                            C_b_o[j,:] = C_b_o[j,:] - 0.5*delta_ij*C_b_o[i,:]
                            C_b_o[i,:] = C_b_o[i,:] - 0.5*delta_ij*C_b_o[j,:]
                            # print(np.dot(C_b_o[i,:], C_b_o[j,:]))


        #ортогонализация
        # for i in range(3):
        #     for j in range(3):
        #         if i!=j:
        #             delta_ij = np.dot(C_b_o[i,:], C_b_o[:,j].T)
        #             if abs(delta_ij)>10**(-6):
        #                 C_b_o[:,j] = C_b_o[:,j] - 0.5*delta_ij*C_b_o[:,i]
        #                 C_b_o[i,:] = C_b_o[i,:] - 0.5*delta_ij*C_b_o[j,:]
        #а сейчас вычисление скооростей

        #нахождение приращения линейных скоростей методом Рунге-Кутта
        delta_alpha_sym = np.array([[0, -delta_alpha_bz[i+4], delta_alpha_by[i+4]], [delta_alpha_bz[i+4], 0, -delta_alpha_bx[i+4]], [-delta_alpha_by[i+4], delta_alpha_bx[i+4], 0]])
        w_b = np.array([[w_bx[i+4]], [w_by[i+4]], [w_bz[i+4]]])
        k_1 = w_b - delta_alpha_sym*delta_W
        k_2 = w_b - delta_alpha_sym*(delta_W+ h_1*k_1/2)
        k_3 = w_b - delta_alpha_sym*(delta_W+ h_1*k_2/2)
        k_4 = w_b - delta_alpha_sym*(delta_W+ h_1*k_3)
        delta_W = delta_W+ (k_1+2*k_2+2*k_3+k_4)/6

        #пересчет линейных скоростей в опорную систему координат

        delta_W_o = np.dot(C_b_o, delta_W)

        V_ox = V_ox+delta_W_o[0,0]+h_4*((U*np.sin(phi)+omega_oz)*V_oy - (U*np.cos(phi)+omega_oy)*V_oz)
        V_oy = V_oy + delta_W_o[1,0] +h_4*(-(U*np.sin(phi)+omega_oz)*V_ox - omega_ox*V_oz)
        V_oz = V_oz + delta_W_o[2,0] +h_4*((U*np.cos(phi)+omega_oy)*V_ox - omega_ox*V_oy - g) #по прикладному алгоритму (блок схема)

        # # по Егорушкину и Салычеву
        # delta_a_up_c = -V_oy**2/(R_phi+Hight) - V_ox**2/(R_lambda+Hight) - 2*V_ox*U*np.cos(phi)
        # V_oz = V_oz + delta_W_o[2,0] - delta_a_up_c -g - c2*(Hight -  Hight_ist)
        # Hight = Hight + (V_oz - c1*(Hight -  Hight_ist))*h_4

        V_x.append(V_ox)
        V_y.append(V_oy)


        omega_ox = -V_oy/R_phi
        omega_oy = U*np.cos(phi)+ V_ox/(R_lambda+Hight)
        omega_oz = U*np.sin(phi)+ V_oy/(R_lambda+Hight)*np.tan(phi)

        #вычисление координат 
        phi = phi + h_6*V_oy/(R_phi+Hight)
        lmbda = lmbda + h_6*V_ox/((R_phi+Hight)*np.cos(phi))

        phi_nav.append(phi*rad_to_deg)
        lmbda_nav.append(lmbda*rad_to_deg)

        psi.append(np.arctan(C_b_o[1,0]/C_b_o[1,1])*rad_to_deg)
        gamma.append(np.arctan(-C_b_o[2,0]/C_b_o[2,2])*rad_to_deg)
        c_0 = np.sqrt(C_b_o[0,1]**2+C_b_o[1,1]**2)
        theta.append(np.arctan(C_b_o[2,1]/c_0)*rad_to_deg)
        flag = True

#построение графиков

# pyplot.plot(V_x)
pyplot.plot(V_y)
pyplot.show()
# pyplot.plot()