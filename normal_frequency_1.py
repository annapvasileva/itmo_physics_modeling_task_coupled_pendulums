from decimal import Decimal, getcontext
import matplotlib.pyplot as plt
from math import pi
from scipy.signal import find_peaks
import numpy as np
from scipy.integrate import solve_ivp

PHI_0 = Decimal(pi / 4)


def sin_decimal(x: Decimal):
    TWO_PI = Decimal(str(2 * pi))
    # Нормализация аргумента к диапазону [-pi, pi]
    x = x % TWO_PI
    if x > Decimal(pi):
        x -= TWO_PI
    
    tmp = getcontext().prec
    getcontext().prec = tmp + 10  # увеличиваем немного, не чрезмерно
    
    n = x
    r = Decimal(0)
    i = 1
    while abs(n) > Decimal("1E-" + str(getcontext().prec)):
        r += n
        n *= Decimal(-1) * x * x / ((2 * i) * (2 * i + 1))
        i += 1
    
    getcontext().prec = tmp
    return r

def f(y1, y2, w0, beta):
    return -w0 * w0 * sin_decimal(y1) - Decimal(2) * beta * y2

def get_next_y2(y1_i : Decimal, y2_i : Decimal, w0 : Decimal, beta : Decimal, h : Decimal):
    k1 = f(y1_i, y2_i, w0, beta)
    k2 = f(y1_i + h / Decimal(2), y2_i + h * k1 / Decimal(2), w0, beta)
    k3 = f(y1_i + h / Decimal(2), y2_i + h * k2 / Decimal(2), w0, beta)
    k4 = f(y1_i + h, y2_i + h * k3, w0, beta)
    
    delta_y2_i = h / Decimal(6) * (k1 + Decimal(2) * k2 + Decimal(2) * k3 + k4)
    
    return y2_i + delta_y2_i
    
def get_w1(natural_frequency : Decimal, beta : Decimal, phi_0 : Decimal = PHI_0, w_0 : Decimal = Decimal(0), t_0 : Decimal = Decimal(0), t_n : Decimal = Decimal(50), t_step : Decimal = Decimal(0.01)):
    steps = int((t_n - t_0) / t_step)
    t_values = [t_step * Decimal(i) for i in range(steps)]
    y1_values = [phi_0]
    y2_values = [w_0]
    
    y1_i = phi_0
    y2_i = w_0
    
    for _ in range(1, steps):
        y2_i_derivative = f(y1_i, y2_i, natural_frequency, beta)
        
        y2_next = get_next_y2(y1_i, y2_i, natural_frequency, beta, t_step)
        
        k1 = y2_i
        k2 = y2_i + t_step * y2_i_derivative / Decimal(2)
        k3 = y2_i + t_step * y2_i_derivative / Decimal(2)
        k4 = y2_i + t_step * y2_i_derivative
        y1_next = y1_i + t_step/Decimal(6)*(k1 + Decimal(2)*k2 + Decimal(2)*k3 + k4)
        
        y1_i, y2_i = y1_next, y2_next
        y1_values.append(y1_i)
        y2_values.append(y2_i)
    
    # -------------------Графики---------------------
    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    # # Первый график: φ(t) [радианы]
    # ax1.plot(t_values, y1_values, label='φ(t)', color='blue')
    # ax1.set_ylabel('φ(t), [рад]', fontsize=12)
    # ax1.set_xlabel('Время, [с]', fontsize=12)
    # ax1.grid(True)
    # ax1.legend()
    # # Второй график: dφ/dt(t) [рад/с]
    # ax2.plot(t_values, y2_values, label='dφ/dt(t)', color='red')
    # ax2.set_ylabel('dφ/dt(t), [рад/с]', fontsize=12)
    # ax2.set_xlabel('Время, [с]', fontsize=12)
    # ax2.grid(True)
    # ax2.legend()
    # plt.tight_layout()
    # plt.show()
    
    # ----------------Поиск частоты--------------------
    t_values_float = [float(t) for t in t_values]
    y1_values_float = [float(y) for y in y1_values]
    # Находим индексы максимумов
    peaks, _ = find_peaks(y1_values_float)
    # Время, соответствующее максимумам
    t_peaks = [t_values_float[i] for i in peaks]
    periods = [t_peaks[i+1] - t_peaks[i] for i in range(len(t_peaks)-1)]
    t_avg = np.mean(periods)
    
    res_w = Decimal(2 * pi) / Decimal(t_avg)
    
    return res_w
    
def get_w1_euler(natural_frequency : Decimal, beta : Decimal, phi_0 : Decimal = PHI_0, w_0 : Decimal = Decimal(0), t_0 : Decimal = Decimal(0), t_n : Decimal = Decimal(50), t_step : Decimal = Decimal(0.01)):
    steps = int((t_n - t_0) / t_step)
    t_values = [t_step * Decimal(i) for i in range(steps)]
    y1_values = [phi_0]
    y2_values = [w_0]
    
    y1_i = phi_0
    y2_i = w_0
    
    for _ in range(1, steps):
        y2_i_derivative = f(y1_i, y2_i, natural_frequency, beta)
        y1_next = y1_i + y2_i * t_step + y2_i_derivative * t_step * t_step / 2
        y2_next = y2_i + y2_i_derivative * t_step
        y1_i = y1_next
        y2_i = y2_next
        y1_values.append(y1_i)
        y2_values.append(y2_i)
    
    # График
    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    # # Первый график: φ(t) [радианы]
    # ax1.plot(t_values, y1_values, label='φ(t)', color='blue')
    # ax1.set_ylabel('φ(t), [рад]', fontsize=12)
    # ax1.set_xlabel('Время, [с]', fontsize=12)
    # ax1.grid(True)
    # ax1.legend()
    # # Второй график: dφ/dt(t) [рад/с]
    # ax2.plot(t_values, y2_values, label='dφ/dt(t)', color='red')
    # ax2.set_ylabel('dφ/dt(t), [рад/с]', fontsize=12)
    # ax2.set_xlabel('Время, [с]', fontsize=12)
    # ax2.grid(True)
    # ax2.legend()
    # plt.tight_layout()
    # plt.show()
    
    t_values_float = [float(t) for t in t_values]
    y1_values_float = [float(y) for y in y1_values]
    # Находим индексы максимумов
    peaks, _ = find_peaks(y1_values_float)
    # Время, соответствующее максимумам
    t_peaks = [t_values_float[i] for i in peaks]
    periods = [t_peaks[i+1] - t_peaks[i] for i in range(len(t_peaks)-1)]
    t_avg = np.mean(periods)
    
    res_w = Decimal(2 * pi) / Decimal(t_avg)
    
    return res_w


# l = Decimal(0.79)
# beta = Decimal(0)
# phi_0 = Decimal(0.17)
# get_w1_euler((Decimal(9.8)/l).sqrt(), beta, phi_0)
