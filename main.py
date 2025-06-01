from decimal import Decimal
from math import pi

from main_graphs import build_graphs, build_graphs_euler
from normal_frequency_1 import get_w1
from normal_frequency_2 import get_w2, G

PHI_0 = Decimal(pi / 4)
DELTA = Decimal(1e-6)

def is_non_negative(x : Decimal):
    return x > -DELTA

def main():
    print("--------Ввод параметров системы--------")
    
    print("L (м): ", end="")
    is_valid = False
    while(not is_valid):
        l = Decimal(input())
        if (l > Decimal(0)):
            is_valid = True
        else:
            print("Число должно быть положительным!")
            
    print("L1 (м): ", end="")
    is_valid = False
    while(not is_valid):
        l1 = Decimal(input())
        if (l1 > Decimal(0) and l1 - l <= DELTA):
            is_valid = True
        else:
            print("Число L1 должно быть положительным и не должно превышать L!")
    
    print("Коэффициент затухания β (с^-1): ", end="")
    is_valid = False
    while(not is_valid):
        beta = Decimal(input())
        if (is_non_negative(beta)):
            is_valid = True
        else:
            print("Число должно быть неотрицательным!")
            
    print("Коэффициент жёсткости пружины k (Н/м): ", end="")
    is_valid = False
    while(not is_valid):
        k = Decimal(input())
        if (is_non_negative(k)):
            is_valid = True
        else:
            print("Число должно быть неотрицательным!")
    print("Начальное отклонение первого маятника φ01 (рад): ", end="")
    is_valid = False
    while(not is_valid):
        phi1_0 = Decimal(input())
        if (-Decimal(pi) - DELTA < phi1_0 < Decimal(pi) + DELTA):
            is_valid = True
        else:
            print("Число должно быть в пределах [-pi/2; pi/2]!")
            
    print("Начальное отклонение второго маятника φ02 (рад): ", end="")
    is_valid = False
    while(not is_valid):
        phi2_0 = Decimal(input())
        if (-Decimal(pi) - DELTA < phi2_0 < Decimal(pi) + DELTA):
            is_valid = True
        else:
            print("Число должно быть в пределах [-pi/2; pi/2]!")
    
    print("Расстояние между точками подвеса маятников d (м): ", end="")
    is_valid = False
    while(not is_valid):
        d = Decimal(input())
        if (d >= Decimal(2) * l):
            is_valid = True
        else:
            print("Число должно быть не меньше 2L!")
    
    print("Масса груза на конце каждого маятника m (кг): ", end="")
    is_valid = False
    while(not is_valid):
        m = Decimal(input())
        if (m > Decimal(0)):
            is_valid = True
        else:
            print("Число должно быть положительным!")
    
    # --------Пример входных данных--------
    # L = 0.79 m
    # L1 = 0.35 m
    # β = 0.146 s^-1
    # k = 10000 N/m
    # φ01 = 0.17 rad
    # φ02 = 0.34 rad
    # d = 2 m
    # m = 50 kg
    
    print("--------Нормальные частоты--------")
    print("ω1 =", get_w1((G / l).sqrt(), beta))
    print("ω2 =", get_w2(k, m, l1, l, beta))
    build_graphs(beta, phi1_0, phi2_0, k, d, l1, m, l)


if __name__ == '__main__':
    main()