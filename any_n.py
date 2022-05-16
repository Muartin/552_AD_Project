from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np
import os

# defining the system of differential equations from 4.1 for any n
def system(t, w, k, kstar, lmbda, M):
        derivatives = []
        n = len(w)      

        #calculate for change in number of monomer w.r.t. time
        dot_w = 0
        for j in range(n - 1):
                dot_w -= k*w[0]*w[j]
        dot_w += -kstar*w[0]*w[-1] - M*w[0] + lmbda
        derivatives.append(dot_w)

        #calculate for change in number of soluble oligomers w.r.t time
        for s in range(1, n - 1):
                dot_w = 0
                for i in range(n):
                        for j in range(n):
                                if (i + 1) + (j + 1) == (s + 1):
                                        dot_w += 0.5*k*w[i]*w[j]
                
                for j in range(n - 1):
                        dot_w -= k*w[s]*w[j]
                
                dot_w += -kstar*w[s]*w[-1] - M*w[s]
                derivatives.append(dot_w)

        #calcualate for change in number of senile plaques w.r.t time
        dot_w = 0
        for i in range(n - 1):
                for j in range(n - 1):
                        if (i + 1) + (j + 1) >= n:
                                dot_w += 0.5*k*w[i]*w[j]
        dot_w -= M*w[-1]
        derivatives.append(dot_w)
        
        return tuple(derivatives)


n = int(input('input a value for n: '))

# the final time and time step for evaluating
tf = 500
dt = 0.1

# initial point to solve with, all initial values are set to 0
initial = tuple([0]*n)

# defining constants
k = 10**(-4)
kstar = 5 * 10**(-6)
lmda = 2
M = 10**(-2)

constants = (k, kstar, lmda, M)

sol = integrate.solve_ivp( system, (0, tf), initial, t_eval=np.arange(0, tf + dt, dt), args=constants)


# plotting as a function of t
plt.title('Plot for n=' + str(n))
for i in range(n):
        this_label = 'Y' + str(i + 1)
        plt.plot(sol.t, sol.y[i], label=this_label)
plt.xlabel('Time (days)')
plt.legend()

filepath = 'figures/'
os.makedirs(os.path.dirname(filepath), exist_ok=True)
plt.savefig(filepath + 'n=' + str(n))
plt.show()