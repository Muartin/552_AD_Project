# imports and such
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


# defining the system of differential equations from 4.1 for n=4
def system(t: int, initial: tuple[int, list, int], k: float, kstar:float , lmda:float, M: float):
        X, Y1, Y2, Z = initial
        dot_X = -k*X**2 - k*X*Y1 -k*X*Y2 - kstar*X*Z - M*X + lmda
        dot_Y1 = 0.5*k*X**2 - k*X*Y1 - k*Y1**2 - k*X*Y2  - kstar*Y1*Z - M*Y1
        dot_Y2 = k*X*Y1 - k*X*Y2 - k*Y1*Y2 - k*Y2**2-k*Y2*Z-M*Y2 
        dot_Z = k*X*Y2 + 0.5*k*Y1**2 + k*Y1*Y2 + 0.5*k*Y2**2 - M*Z
        return dot_X, dot_Y1, dot_Y2, dot_Z


# the final time and time step for evaluating
tf = 500
dt = 0.1

# initial point to solve with
initial = (0, 0, 0, 0)

# defining constants
k = 10**(-4)
kstar = 5 * 10**(-6)
lmda = 2
M = 10**(-2)

constants = (k, kstar, lmda, M)

sol = integrate.solve_ivp( system, (0, tf), initial, t_eval=np.arange(0, tf + dt, dt), args=constants )


# plotting as a function of t
plt.plot(sol.t, sol.y[0], label='X')
plt.plot(sol.t, sol.y[1], label='Y1')
plt.plot(sol.t, sol.y[2], label='Y2')
plt.plot(sol.t, sol.y[3], label='Z')
plt.xlabel('Time (days)')
plt.legend()
plt.show()