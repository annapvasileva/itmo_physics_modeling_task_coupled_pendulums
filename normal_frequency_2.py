from decimal import Decimal, getcontext
import matplotlib.pyplot as plt
from math import pi
from scipy.signal import find_peaks
import numpy as np
from normal_frequency_1 import sin_decimal, PHI_0

G = Decimal(9.8)            # m/s^2

def f(y1 : Decimal, y2 : Decimal, beta : Decimal, k : Decimal, m : Decimal, l1 : Decimal, l : Decimal):
    w0 = (G / l).sqrt()
    
    return -w0 * w0 * sin_decimal(y1) - Decimal(2) * beta * y2 - k / m * (l1 * l1 / l / l) * sin_decimal(Decimal(2) * y1)

def get_next_y2(y1_i : Decimal, y2_i : Decimal, beta : Decimal, h : Decimal, k : Decimal, m : Decimal, l1 : Decimal, l : Decimal):
    k1 = f(y1_i, y2_i, beta, k, m, l1, l)
    k2 = f(y1_i + h / Decimal(2), y2_i + h * k1 / Decimal(2), beta, k, m, l1, l)
    k3 = f(y1_i + h / Decimal(2), y2_i + h * k2 / Decimal(2), beta, k, m, l1, l)
    k4 = f(y1_i + h, y2_i + h * k3, beta, k, m, l1, l)
    
    delta_y2_i = h / Decimal(6) * (k1 + Decimal(2) * k2 + Decimal(2) * k3 + k4)
    
    return y2_i + delta_y2_i
    
def get_w2(k : Decimal, m : Decimal, l1 : Decimal, l : Decimal, beta : Decimal, phi_0 : Decimal = PHI_0, w_0 : Decimal = Decimal(0), t_0 : Decimal = Decimal(0), t_n : Decimal = Decimal(50), t_step : Decimal = Decimal(0.01)):
    steps = int((t_n - t_0) / t_step)
    t_values = [t_step * Decimal(i) for i in range(steps)]
    y1_values = [phi_0]
    y2_values = [w_0]
    
    y1_i = phi_0
    y2_i = w_0
    
    for _ in range(1, steps):
        y2_i_derivative = f(y1_i, y2_i, beta, k, m, l1, l)
        
        y2_next = get_next_y2(y1_i, y2_i, beta, t_step, k, m, l1, l)
        
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
    
def get_w2_euler(k : Decimal, m : Decimal, l1 : Decimal, l : Decimal, beta : Decimal, phi_0 : Decimal = PHI_0, w_0 : Decimal = Decimal(0), t_0 : Decimal = Decimal(0), t_n : Decimal = Decimal(50), t_step : Decimal = Decimal(0.01)):
    steps = int((t_n - t_0) / t_step)
    t_values = [t_step * Decimal(i) for i in range(steps)]
    y1_values = [phi_0]
    y2_values = [w_0]
    
    y1_i = phi_0
    y2_i = w_0
    
    for _ in range(1, steps):
        y2_i_derivative = f(y1_i, y2_i, beta, k, m, l1, l)
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
    
