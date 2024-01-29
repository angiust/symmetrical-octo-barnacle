import numpy as np
import sys

def my_print(x, v, t):
    print(f"{x[0]} {x[1]} {x[2]} {v[0]} {v[1]} {v[2]} {t}") # così mi stampa senza parentesi quadre che danno fastidio a gnuplot

def euler(x, v, t, gradient, dt):
    return x + dt * v, v - dt * gradient(x, v, t)

def loop(x, v, gradient, dt, t):  # seguo il moto da 0 a t
    n = 0
    while dt * n <= t:
        if n % 1_000 == 0: # ogni quanto voglio stampare
            my_print(x, v, dt * n)
        x, v = euler(x, v, dt * n, gradient, dt)
        n += 1

if __name__ == "__main__": # parse degli argomenti da terminale
    if len(sys.argv) != 9: # voglio 9 elementi: nome_del_programma.py x0 y0 z0 vx0 vy0 vz0 dt t
        print("Usage: python3 nome_del_programma.py x0 y0 z0 vx0 vy0 vz0 dt t")
        sys.exit(1)

    x_0 = np.array([float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])])
    v_0 = np.array([float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])])
    dt = float(sys.argv[7])
    total_time = float(sys.argv[8])
# parametri del modello e dati iniziali
a = -1
b = 0.25
g = 2.5
d = 0.1
o = 2
gradient = lambda x, v, t: d * v + a * x + b * x**3 - g * np.cos(o * t)  # duffing oscillator
loop(x_0, v_0, gradient, dt, total_time) # lancio l'integrazione numerica

