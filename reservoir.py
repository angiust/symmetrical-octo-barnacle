import numpy as np


class Reservoir:
    def __init__(self, n_inputs=1, n_neurons, inp_scaling=1., parameters, delay):
        self.n_inputs = n_inputs
        self.n_neurons = n_neurons
        self.inp_scaling = inp_scaling
        self.parameters = parameters
        self.delay = delay
        self.leak_rate = np.exp(-np.log(1+delay/n_neurons))

        # initialize mask
        self.mask = np.random.uniform(low=-1., high=1., size=(n_neurons, n_inputs)) * inp_scaling

        # initialize previous states
        self.prev_states = np.zeros((n_neurons,))
        #initialize states
        self.states = []

        def kernel_function(self, x_prev, u):
        """
        kernel function of the reservoir, its parameters:
        x_prev: previous state
        it's' possible to add x: current state
        u: current input
        to add theta: parameters of the kernel (eg: in the case of Mackey-Glass kernel (gamma, eta, p))
        """
        # Mackey-Glass non linear kernel:  eta*(x+gamma*u)/(1+(x+gamma*u)**p)
        # now is implemented the hyperbolic tangent
        return np.tanh(x_prev + u)

    def forward(self, input_signal, wout):
        """
        execution of a temporal step of the reservoir: returns current states of neurons and the output, updates previous state
        u: input signal
        wout: output weights
        """
        # apply the mask to the input
        input_signal *= self.mask
        # initialize state of neurons
        x = np.zeros((self.n_neurons,))
        # loop for every neurons
        for i in range(self.n_neurons):
            if i == 0: # the first neuron depends on its values in the previous state and on the values of the last neuron in the previous state
                x[i] = self.leak_rate*self.prev_states[self.n_neurons-1] + (1-self.leak_rate)*self.kernel_function(self.prev_states[i], input_signal[i])
            else: # others neurons depend on its values in the previous state and on the values of the previous neuron in the current state
                x[i] = self.leak_rate*x[i-1] + (1-self.leak_rate)*self.kernel_function(self.prev_states[i], input_signal[i])
        """ it is possible do it without saving previous states and using self.states:

        for i in range(self.n_neurons):
            if i == 0:
                x[i] = self.leak_rate * self.states[-1][i] + (1 - self.leak_rate) * self.kernel_function(self.states[-1][i], input_signal[i])
            else:
                x[i] = self.leak_rate * x[i-1] + (1 - self.leak_rate) * self.kernel_function(self.states[-1][i], input_signal[i])
        """
        # update previous state
        self.prev_states = x.copy()
        self.states.append(x.copy())
        # compute the output
        y = x @ wout

        return x, y

#to do: the training of wout, this is why self.states is for
