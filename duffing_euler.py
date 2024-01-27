import numpy as np

def my_print(x, v, t):
    print(f"{x} {v} {t}")

def euler(x, v, t, gradient, dt):
    return x + dt * v, v - dt * gradient(x, v, t)

def loop(x, v, gradient, dt, t): # seguo il moto da 0 a t
    n = 0
    while(dt * n <= t):
        if n % 1_000 == 0: # ogni quanto voglio stampare
            my_print(x, v, dt * n)
        x, v = euler(x, v, dt * n, gradient, dt)
        n += 1
#parametri del modello e dati iniziali
a = -1
b = 0.25
g = 2.5;
d = 0.1
o = 2
x_0 = np.array([0,0,0])
v_0 = np.array([1,0,0])
dt = 0.000001
gradient = lambda x,v,t : d*v + a*x + b*x**3 - g* np.cos(o*t) # duffing oscillator
loop(x_0, v_0, gradient, dt, 100)
