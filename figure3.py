# imports and such
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


# defining the system of differential equations from 4.1 for n=4
def system(t: int, initial: tuple[int, list, int], k: float, kstar:float , lmda:float, M: float):
        Y1, Y2, Y3, Y4, Y5 = initial
        dot_Y1 = -k*Y1**2 - k*Y1*Y2-k*Y1*Y3 - k*Y1*Y4 - kstar*Y1*Y4 - M*Y1 + lmda
        dot_Y2 = 0.5*k*Y1**2 - k*Y1*Y2 - k*Y2**2 - k*Y2*Y3 - k*Y2*Y4 - kstar*Y2*Y5 - M*Y2
        dot_Y3 = k*Y1*Y2 - k*Y1*Y3 - k*Y2*Y3 - k*Y3**2 - k*Y3*Y4 - kstar*Y3*Y5 - M*Y3
        dot_Y4 = k*Y1*Y3 + 0.5*k*Y2**2 - k*Y1*Y4 - k*Y2*Y4 - k*Y3*Y4 - k*Y4**2 - kstar*Y4*Y5 - M*Y4
        dot_Y5 = k*Y1*Y4 + k*Y2*Y3 + 0.5*k*Y3**2 + k*Y3*Y4 + 0.5*k*Y4**2 - M*Y5 

        return dot_Y1, dot_Y2, dot_Y3, dot_Y4, dot_Y5


# the final time and time step for evaluating
tf = 500
dt = 0.1

# initial point to solve with
initial = (0, 0, 0, 0, 0)

# defining constants
k = 10**(-4)
kstar = 5 * 10**(-6)
lmda = 2
M = 10**(-2)

constants = (k, kstar, lmda, M)

sol = integrate.solve_ivp( system, (0, tf), initial, t_eval=np.arange(0, tf + dt, dt), args=constants )


# plotting as a function of t
plt.plot(sol.t, sol.y[0], label='Y1')
plt.plot(sol.t, sol.y[1], label='Y2')
plt.plot(sol.t, sol.y[2], label='Y3')
plt.plot(sol.t, sol.y[3], label='Y4')
plt.plot(sol.t, sol.y[4], label='Y5')
plt.xlabel('Time (days)')
plt.legend()
plt.show()