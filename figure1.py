# imports and such
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


# defining the system of differential equations from 4.1 for n=3
def system(t, initial, k, kstar, lmda, M1, M2, M3):
        X, Y, Z = initial
        dot_X = -k*X**2 - k*X*Y - kstar*X*Z - M1*X + lmda
        dot_Y = 0.5*k*X**2 - k*X*Y - k*Y**2 - kstar*Y*Z - M2*Y
        dot_Z = 0.5*k*Y**2 + k*X*Y - M3*Z
        return dot_X, dot_Y, dot_Z


# the final time and time step for evaluating
tf = 500
dt = 0.1

# initial point to solve with
initial = (0, 0, 0)

# defining constants
k = 10**(-4)
kstar = 5 * 10**(-6)
lmda = 2
M1 = 10**(-2)
M2 = 10**(-2)
M3 = 10**(-2)
constants = (k, kstar, lmda, M1, M2, M3)

sol = integrate.solve_ivp( system, (0, tf), initial, t_eval=np.arange(0, tf + dt, dt), args=constants )


# plotting as a function of t
plt.plot(sol.t, sol.y[0], label='X')
plt.plot(sol.t, sol.y[1], label='Y')
plt.plot(sol.t, sol.y[2], label='Z')
plt.xlabel('Time (days)')
plt.legend()
plt.show()