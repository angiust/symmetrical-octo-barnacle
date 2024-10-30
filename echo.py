import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # o 'Qt5Agg' se hai PyQt5 installato
import matplotlib.pyplot as plt

# it's possible to make the random numbers reproducible by setting the seed 
# np.random.seed(42)

# parameters
input_size = 20 # we will refer to it as N_u for the sake of brevity
target_size = 3 # N_y
reservoir_size = 100 # N_x
spectral_radius = 0.9 # rho , as a guiding principle, ρ(W) should be set greater for tasks where a more extensive history of the input is required to perform it, and smaller for tasks where the current output y(n) depends more on the recent history of u(n). [Lukoševičius, 2012]
leaking_rate = 0.5 # also called alpha in some references, it could also be adapted online
beta = 0.1 # it's parameter for the ridge regression, it balances the trade-off between minimizing the prediction error and controlling the norm of the output weights

#initialize input and target vectors, this is just for now, we will have to put some sensible values
u = np.random.rand(input_size) # it is important to have a good input scaling, as it will affect the learning process
target = np.random.rand(target_size)

# initialize status
x = np.zeros((reservoir_size, 10)) 
# np.column_stack((matrice, vettore_uno_dx))

# initialize input weights, a matrix of dimension N_x x (1+N_u)
W_in = np.random.rand(reservoir_size, 1 + input_size)
# initialize reservoir weights, a matrix of dimension N_x x (1+N_u)
W = np.random.rand(reservoir_size, reservoir_size)
# adjust W in order to have a good spectral radius for the esp
radius = np.max(np.abs(np.linalg.eigvals(W)))
W = W * spectral_radius / radius

# initialize output weights, a matrix of dimension N_y x (1+N_u+N_x)
W_out = np.random.rand(target_size, 1 + input_size + reservoir_size)

for i in range(10):
    # u is the input vector
    u_new = np.random.rand(input_size)
    #print(u_new)
    u = np.column_stack((u, u_new))
    # target is the target vector
    target_new = np.random.rand(target_size)
    target = np.column_stack((target, target_new))

    # update the reservoir state
    x[:,i] = (1 - leaking_rate) * x[:,i-1] + leaking_rate * np.tanh(np.dot(W, x[:,i-1]) + np.dot(W_in, np.concatenate(([1], u[:,i]))))
    #print(f"x[:,{i}]:")
    #print(x[:,i])

    # compute the output
    y = np.dot(W_out, np.concatenate((np.array([1]), u[:,i], x[:,i])))
    #print("y:")
    #print(y)

    # update the output weights
    #W_out = W_out - np.outer(y - target[:,i], np.concatenate((np.array([1]), u[:,i], x[:,i])))

    # update the reservoir weights
    #W = W - np.outer(x, np.dot(W_out[:, 1:], y - target[:,i]))

    # plot the output
    plt.plot(y[0])
    plt.plot(target[0])
    plt.show()
    plt.savefig("nomefile.png")


