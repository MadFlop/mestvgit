import numpy as np
import matplotlib.pyplot as plt


# Прочностные характеристики
B = 25.0            # Класс бетона по прочности на сжатие в МПа
Rb = -36.0          # Предел прочности бетона при сжатии в МПа
Rbt = 3.0           # Предел прочности бетона при растяжении в МПа
Rs = 400.0          # Предел прочности арматуры при растяжении и сжатии в МПа
Rb_ser = Rb
Rbt_ser = Rbt


# Жесткостные характеристики
Eb = 3e4            # Начальный модуль упругости бетона в МПа
Es = 2e5            # Модуль упругости стали в МПа
E_tau_s = 1500.0    # 
eps_y = Rs/Es       # Предел пропорциональности стали


# Геометрические характеристики (ед. изм - метр)
b = 0.1
h = 0.15
a = 0.02


# Количество расчетных слоев бетона и стержней арматуры
h_b = 0.003         # в метрах
n_b = h/h_b
ib = np.arange(n_b)
jb = np.array([0,1])
print(n_b,ib,jb)


temp = 20           # Температура бетона
lambda_ = 1.0       # Коэффициент для определения eps'_b, зависящий от вида бетона


# Деформации в вершине диаграммы работы бетона при temp = +20
eps1_b = -B/Eb * lambda_ * (1 + (0.8-0.15*B**2/1e4) * lambda_*B/60+0.2*lambda_/B) / (0.12+1.03*B/60+0.2/B)


#При определении коэффициента Мурашева psi_sj при кратковременном действии нагрузки
fi_sl = 1.0


#Функция диаграммы Карпенко от eps_b
def karpenkoeps(eps_b):
    if eps_b == 0:
        sigma_b = 0.0
    else:
        beta_temp_E = 1 + 0.2 * (20-temp)/90
        beta_temp_eps = 1 + 0.55 * (20-temp)/90
        if eps_b < 0:
            beta_temp_R = 1 + 0.6 * (20-temp)/90
            sigma1_b = Rb_ser
            sigma1_btemp = sigma1_b * beta_temp_R
            nu1_b = sigma1_btemp / (eps1_b*beta_temp_eps*Eb*beta_temp_E)
            eps1_btemp = eps1_b * beta_temp_eps
        elif eps_b > 0:
            beta_temp_Rt = 1 + 1.3*(20-temp)/90
            sigma1_b = Rbt_ser
            sigma1_btemp = sigma1_b*beta_temp_Rt
            nu1_b = 0.6 + 0.15 * sigma1_btemp/2.5
            eps1_btemp = sigma1_btemp / (nu1_b*beta_temp_eps*Eb*beta_temp_E)
        if 0 < abs(eps_b) < abs(eps1_btemp):
            nu0 = 1.0
            omega1 = 2 - 2.5*nu1_b
        elif abs(eps_b) > abs(eps1_btemp):
            nu0 = 2.05 * nu1_b
            omega1 = 1.95 * nu1_b - 0.138
        u = eps_b*Eb*beta_temp_E/sigma1_btemp
        nu_b = (-(((nu0**4-4*nu1_b*nu0**3+6*nu1_b**2*nu0**2-4*nu1_b**3*nu0+nu1_b**4)*omega1**2 + 
        (16*nu1_b*nu0**3-4*nu0**4-20*nu1_b**2*nu0**2+8*nu1_b**3*nu0)*omega1 + 
        (4*nu0**4-16*nu1_b*nu0**3+20*nu1_b**2*nu0**2-8*nu1_b**3*nu0))*u**2 + 
        (8*nu1_b**2*nu0-4*nu1_b*nu0**2-4*nu1_b**3)*omega1*u+4*nu0**2-8*nu1_b*nu0+4*nu1_b**2)**0.5 + 
        ((nu0**2-2*nu1_b*nu0+nu1_b**2)*omega1*u-2*nu1_b)) / (((2*nu0**2-4*nu1_b*nu0+2*nu1_b**2)*omega1 + (4*nu1_b*nu0-2*nu0**2-2*nu1_b**2))*u**2 - 2.0)
        sigma_b = eps_b*Eb*beta_temp_E*nu_b
    return sigma_b


eps_b1 = np.array([-0.005111863925,-0.004039917307,-0.003421459796,
-0.002995680928,-0.002573376695,-0.002227891366,-0.001810110132,
-0.001466967584,-0.001170221538,-0.0009070284521,-0.0007056858944,
-0.0004698355557,-0.0002835760761,-0.0001300544522,0.0])
sigma_b1 = np.zeros(eps_b1.size)
eps_b2 = np.array([0.0,2.843*1e-5,5.726*1e-5,8.696*1e-5,1.182*1e-4,1.404*1e-4,1.644*1e-4,
1.859*1e-4,2.047*1e-4,2.206*1e-4,2.294*1e-4,2.363*1e-4,2.543*1e-4,2.743*1e-4,2.977*1e-4])
sigma_b2 = np.zeros(eps_b2.size)
for i in range(eps_b1.size):
    sigma_b1[i] = karpenkoeps(eps_b1[i])
for i in range(eps_b2.size):
    sigma_b2[i] = karpenkoeps(eps_b2[i])