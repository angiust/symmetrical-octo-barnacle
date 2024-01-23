def sum(a,b):
    return [a[0] + b[0], a[1] + b[1], a[2]+b[2]]

def prod(s,a):
    return [s*a[0], s*a[1], s*a[2]]

def my_print(x, v, t):
    print(f"{x[0]} {x[1]} {x[2]} {v[0]} {v[1]} {v[2]} {t}")

def index(x, dx):
    return [round(x[0]/dx), round(x[1]/dx), round(x[2]/dx)]

def discretize(x, dx):
    return prod(dx, index(x, dx))

def central_dir_deriv(x, f, dx, d):
    return (f(sum(x,prod(dx,d)))-f(sum(x,prod(-dx,d))))/(2*dx)

def fin_diff_grad(x, f, dx):
    return discretize([central_dir_deriv(x,f, dx, [1,0,0]),
        central_dir_deriv(x,f, dx, [0,1,0]),
        central_dir_deriv(x,f, dx, [0,0,1])], dx)

def euler(x, v, f, dx, dt):
    return discretize(sum(x, prod(dt, v)), dx), sum(v, prod(-dt, fin_diff_grad(x, f, dx)))

def loop(x, v, f, dx, dt, t):
    n = 0
    while(dt * n <= t):
        if n % 1_0 == 0:  #il numero per cui e' valutato il resto di n corrisponde a ogni quanti step temporali stampare i valori di posizione, velocita' e tempo
            my_print(x, v, dt * n)
        x, v = euler(x, v, f, dx, dt)
        n += 1


x_0 = [1,0,0]
v_0 = [0,0,0]
dx = 0.0000001
dt = 0.0001
f = lambda x: x[0]*x[0]+x[1]*x[1]+x[2]*x[2]

loop(x_0, v_0, f, dx, dt, dt*1000000)
