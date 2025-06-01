from decimal import Decimal
import matplotlib.pyplot as plt

from normal_frequency_1 import sin_decimal
from normal_frequency_2 import G

W1_0 = Decimal(0)           # s^-1
W2_0 = Decimal(0)           # s^-1

EPSILON = Decimal(0e-7)

def cos_decimal(x: Decimal):
    sin_dec = sin_decimal(x)
    sin_sq = sin_dec * sin_dec
    if sin_sq > Decimal(1):
        sin_sq = Decimal(1)
    return (Decimal(1) - sin_sq).sqrt()

def delta_l(d: Decimal, l1: Decimal, y1: Decimal, y2: Decimal):
    return ((d + l1 * (sin_decimal(y2) - sin_decimal(y1))) ** 2 + l1 * l1 * (cos_decimal(y2) - cos_decimal(y1)) ** 2).sqrt() - d

def o1h(l1: Decimal, d: Decimal, y1: Decimal, y2: Decimal):
    b = l1 * (cos_decimal(y2) - cos_decimal(y1)) / (d + l1 * (sin_decimal(y2) - sin_decimal(y1)))
    j = l1 * l1 * sin_decimal(y1) * (cos_decimal(y2) - cos_decimal(y1)) / (d + l1 * (sin_decimal(y2) - sin_decimal(y1))) + l1 * cos_decimal(y1)
    x_a = l1 * sin_decimal(y1)
    y_a = l1 * cos_decimal(y1)
    x_b = d + l1 * sin_decimal(y2)
    y_b = l1 * cos_decimal(y2)
    x_h = -j * (y_b - y_a) / (x_b - x_a + b * (y_b - y_a))
    y_h = b * x_h + j
    return (x_h * x_h + y_h * y_h).sqrt()

def o2k(l1: Decimal, d: Decimal, y1: Decimal, y2: Decimal):
    b = l1 * (cos_decimal(y2) - cos_decimal(y1)) / (d + l1 * (sin_decimal(y2) - sin_decimal(y1)))
    j = l1 * l1 * sin_decimal(y1) * (cos_decimal(y2) - cos_decimal(y1)) / (d + l1 * (sin_decimal(y2) - sin_decimal(y1))) + l1 * cos_decimal(y1)
    x_a = l1 * sin_decimal(y1)
    y_a = l1 * cos_decimal(y1)
    x_b = d + l1 * sin_decimal(y2)
    y_b = l1 * cos_decimal(y2)
    x_k = (d * (x_b - x_a) - j * (y_b - y_a))/(x_b - x_a + b * (y_b - y_a))
    y_k = b * x_k + j
    return ((x_k - d) * (x_k - d) + y_k * y_k).sqrt()

def f1(y1: Decimal, y2: Decimal, y3 : Decimal, w0: Decimal, beta: Decimal, k: Decimal, d: Decimal, l1: Decimal, m: Decimal, l: Decimal):
    return -w0 * w0 * sin_decimal(y1) - Decimal(2) * beta * y3 + k * delta_l(d, l1, y1, y2) * o1h(l1, d, y1, y2) / m / l / l

def f2(y1: Decimal, y2: Decimal, y4 : Decimal, w0: Decimal, beta: Decimal, k: Decimal, d: Decimal, l1: Decimal, m: Decimal, l: Decimal):
    return -w0 * w0 * sin_decimal(y2) - Decimal(2) * beta * y4 - k * delta_l(d, l1, y1, y2) * o2k(l1, d, y1, y2) / m / l / l

def energy(m : Decimal, l : Decimal, phi : Decimal, phi_derivative : Decimal, k :Decimal, dl : Decimal):
    return m * G * l * (Decimal(1) - cos_decimal(phi)) + m * l * l * phi_derivative * phi_derivative / Decimal(2) + k / Decimal(4) * dl * dl

def energy_too_huge(m : Decimal, l : Decimal, k :Decimal, dl : Decimal, phi : Decimal, phi_derivative : Decimal, delta_energy : Decimal, energy_0 : Decimal):
    current_energy = energy(m, l, phi, phi_derivative, k, dl)
    
    return current_energy > energy_0 + delta_energy + EPSILON

def correct_w(energy_0 : Decimal, work_of_resistance : Decimal, m : Decimal, l : Decimal, phi : Decimal, k : Decimal, dl : Decimal):
    energy = energy_0 + work_of_resistance
    new_phi_derivative = Decimal (2) / m / l / l * energy - Decimal(2) * G / l * (Decimal(1) - cos_decimal(phi)) - k * dl * dl / Decimal(2) / m / l / l
    if new_phi_derivative < EPSILON:
        new_phi_derivative = Decimal(0)
    
    
    return new_phi_derivative.sqrt()

def get_next(y1_i : Decimal, y2_i : Decimal, y3_i : Decimal, y4_i : Decimal, w0 : Decimal, beta : Decimal, h : Decimal, k : Decimal, d : Decimal, l1 : Decimal, m : Decimal, l : Decimal):
    k1_1 = f1(y1_i, y2_i, y3_i, w0, beta, k, d, l1, m, l)
    k2_1 = f2(y1_i, y2_i, y4_i, w0, beta, k, d, l1, m, l)
    k1_2 = f1(y1_i + k1_1 * h / Decimal(2), y2_i + k2_1 * h / Decimal(2), y3_i, w0, beta, k, d, l1, m, l)
    k2_2 = f2(y1_i + k1_1 * h / Decimal(2), y2_i + k2_1 * h / Decimal(2), y4_i, w0, beta, k, d, l1, m, l)
    k1_3 = f1(y1_i + k1_2 * h / Decimal(2), y2_i + k2_2 * h / Decimal(2), y3_i, w0, beta, k, d, l1, m, l)
    k2_3 = f2(y1_i + k1_2 * h / Decimal(2), y2_i + k2_2 * h / Decimal(2), y4_i, w0, beta, k, d, l1, m, l)
    k1_4 = f1(y1_i + k1_3 * h, y2_i + k2_3 * h, y3_i, w0, beta, k, d, l1, m, l)
    k2_4 = f2(y1_i + k1_3 * h, y2_i + k2_3 * h, y4_i, w0, beta, k, d, l1, m, l)

    delta_y3_i = h / Decimal(6) * (k1_1 + Decimal(2) * k1_2 + Decimal(2) * k1_3 + k1_4)
    delta_y4_i = h / Decimal(6) * (k2_1 + Decimal(2) * k2_2 + Decimal(2) * k2_3 + k2_4)
    
    k1 = y3_i
    k2 = y3_i + h * k1_1 / Decimal(2)
    k3 = y3_i + h * k1_1 / Decimal(2)
    k4 = y3_i + h * k1_1
    y1_next = y1_i + h/Decimal(6)*(k1 + Decimal(2)*k2 + Decimal(2)*k3 + k4)
    
    k1 = y4_i
    k2 = y4_i + h * k2_1 / Decimal(2)
    k3 = y4_i + h * k2_1 / Decimal(2)
    k4 = y4_i + h * k2_1
    y2_next = y2_i + h/Decimal(6)*(k1 + Decimal(2)*k2 + Decimal(2)*k3 + k4)
    
    return y1_next, y2_next, y3_i + delta_y3_i, y4_i + delta_y4_i
    
def build_graphs(beta : Decimal, phi1_0 : Decimal, phi2_0 : Decimal, k : Decimal, d : Decimal, l1 : Decimal, m : Decimal, l : Decimal, t_0 : Decimal = Decimal(0), t_n : Decimal = Decimal(50), t_step : Decimal = Decimal(0.01)):
    natural_frequency = (G / l).sqrt()
    
    steps = int((t_n - t_0) / t_step)
    t_values = [t_step * Decimal(i) for i in range(steps)]
    y1_values = [phi1_0]
    y2_values = [phi2_0]
    y3_values = [W1_0]
    y4_values = [W2_0]
    
    delta_l_0 = delta_l(d, l1, phi1_0, phi2_0)
    energy_0_1 = m * G * l * (Decimal(1) - cos_decimal(phi1_0)) + k / 4 * delta_l_0 * delta_l_0
    energy_0_2 = m * G * l * (Decimal(1) - cos_decimal(phi2_0)) + k / 4 * delta_l_0 * delta_l_0
    c = Decimal(2) * m * beta
    work_of_resistance_1 = Decimal(0)
    work_of_resistance_2 = Decimal(0)
    
    y1_i = phi1_0
    y2_i = phi2_0
    y3_i = y3_values[0]
    y4_i = y4_values[0]
    
    for _ in range(1, steps):
        y1_i, y2_i, y3_i, y4_i = get_next(y1_i, y2_i, y3_i, y4_i, natural_frequency, beta, t_step, k, d, l1, m, l)

        abs_y3_i = abs(y3_i)
        abs_y4_i = abs(y4_i)
        work_of_resistance_1 -= c * l * l * abs_y3_i * abs_y3_i * t_step
        work_of_resistance_2 -= c * l * l * abs_y4_i * abs_y4_i * t_step
        delta_l_i = delta_l(d, l1, y1_i, y2_i)
        if energy_too_huge(m, l, k, delta_l_i, y1_i, y3_i, work_of_resistance_1, energy_0_1):
            if y3_i < EPSILON:
                y3_i = -correct_w(energy_0_1, work_of_resistance_1, m, l, y1_i, k, delta_l_i)
            else:
                y3_i = correct_w(energy_0_1, work_of_resistance_1, m, l, y1_i, k, delta_l_i)
                
            
        if energy_too_huge(m, l, k, delta_l_i, y2_i, y4_i, work_of_resistance_2, energy_0_2):
            if y4_i < EPSILON:
                y4_i = -correct_w(energy_0_2, work_of_resistance_2, m, l, y2_i, k, delta_l_i)
            else:
                y4_i = correct_w(energy_0_2, work_of_resistance_2, m, l, y2_i, k, delta_l_i)
                
            
        y1_values.append(y1_i)
        y2_values.append(y2_i)
        y3_values.append(y3_i)
        y4_values.append(y4_i)
    
    # -------------------Графики---------------------
    _, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    # Графики: φ1(t) и φ2(t)  [радианы]
    ax1.plot(t_values, y1_values, label='φ1(t)', color='cyan')
    ax1.plot(t_values, y2_values, label='φ2(t)', color='purple')
    ax1.set_ylabel('φ(t), [рад]', fontsize=12)
    ax1.set_xlabel('Время, [с]', fontsize=12)
    ax1.grid(True)
    ax1.legend()
    # Графики w1(t) и w2(t) [рад/с]
    ax2.plot(t_values, y3_values, label='ω1(t)', color='cyan')
    ax2.plot(t_values, y4_values, label='ω2(t)', color='purple')
    ax2.set_ylabel('ω(t), [рад/с]', fontsize=12)
    ax2.set_xlabel('Время, [с]', fontsize=12)
    ax2.grid(True)
    ax2.legend()
    plt.tight_layout()
    plt.show()

def build_graphs_euler(beta: Decimal, phi1_0: Decimal, phi2_0: Decimal, k: Decimal, d: Decimal, l1: Decimal, m: Decimal, l: Decimal, t_0: Decimal = Decimal(0), t_n: Decimal = Decimal(50), t_step: Decimal = Decimal(0.01)):
    natural_frequency = (G / l).sqrt()
    
    steps = int((t_n - t_0) / t_step)
    t_values = [t_step * Decimal(i) for i in range(steps)]
    y1_values = [phi1_0]
    y2_values = [phi2_0]
    y3_values = [W1_0]
    y4_values = [W2_0]
    
    y1_i = y1_values[0]
    y2_i = y2_values[0]
    y3_i = y3_values[0]
    y4_i = y4_values[0]
    
    for _ in range(1, steps):
        y3_i_derivative = f1(y1_i, y2_i, y3_i, natural_frequency, beta, k, d, l1, m, l)
        y4_i_derivative = f2(y1_i, y2_i, y4_i, natural_frequency, beta, k, d, l1, m, l)
        y1_next = y1_i + y3_i * t_step + y3_i_derivative * t_step * t_step / Decimal(2)
        y2_next = y2_i + y4_i * t_step + y4_i_derivative * t_step * t_step / Decimal(2)
        y3_next = y3_i + y3_i_derivative * t_step
        y4_next = y4_i + y4_i_derivative * t_step
        y1_i = y1_next
        y2_i = y2_next
        y3_i = y3_next
        y4_i = y4_next
        y1_values.append(y1_i)
        y2_values.append(y2_i)
        y3_values.append(y3_i)
        y4_values.append(y4_i)
    
    # -------------------Графики---------------------
    _, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    # Графики: φ1(t) и φ2(t)  [радианы]
    ax1.plot(t_values, y1_values, label='φ1(t)', color='cyan')
    ax1.plot(t_values, y2_values, label='φ2(t)', color='purple')
    ax1.set_ylabel('φ(t), [рад]', fontsize=12)
    ax1.set_xlabel('Время, [с]', fontsize=12)
    ax1.grid(True)
    ax1.legend()
    # Графики w1(t) и w2(t) [рад/с]
    ax2.plot(t_values, y3_values, label='ω1(t)', color='cyan')
    ax2.plot(t_values, y4_values, label='ω2(t)', color='purple')
    ax2.set_ylabel('ω(t), [рад/с]', fontsize=12)
    ax2.set_xlabel('Время, [с]', fontsize=12)
    ax2.grid(True)
    ax2.legend()
    plt.tight_layout()
    plt.show()
