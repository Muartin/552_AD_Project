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


def degradation(gammas, odeys):
        deg_vals = []
        for t in range(len(odeys[0])):
                deg = 0
                for s in range(1, len(odeys) - 1):
                        print(s)
                        deg += gammas[s]*odeys[s][t]
                deg_vals.append(deg)
        return deg_vals


def ad_progression(ad, tf, deg, deglimit, theta):
        # updated twice a year
        t = 0.5
        while t <= tf:
                deg_val = int(t*400)
                if deg[deg_val] > deglimit:
                        new_ad = ad[len(ad - 1)] + theta * (deg[deg_val] - dlimit)
                        ad.append(new_ad)
                else:
                        ad.append(ad[len(ad)-1])
                print("ad(" + str(t) + ") = " + str(ad[len(ad)-1]))
                t += 0.5
        return ad


n = int(input('input a value for n: '))

# the final time and time step for evaluating
tf = 20
dt = 0.05

# initial point to solve with, all initial values are set to 0
initial = tuple([0]*n)

# defining constants
k = 10**(-4)
kstar = 5 * 10**(-6)
lmda = 2.548
M = 10**(-2)

constants = (k, kstar, lmda, M)

sol = integrate.solve_ivp( system, (0, tf), initial, t_eval=np.arange(0, tf + dt, dt), args=constants)

# sol.y returns a list of matrices of y values for each y
# sol.t returns a list of time points of t

# set gamma influences to be uniform among oligomers
gmma = [0]
for i in range(n-2):
        gmma.append(1/(n-2))
gmma.append(0)
print("gammas: " + str(gmma))
# calculate degradation
d = degradation(gmma, sol.y)
max = 0
# grab maximum value of degradation
for degrade in d:
        if degrade > max:
                max = degrade
print("max deg = " + str(max))

# calculate ad progression
dlimit = 22
thta = 10**(-3)
ad = [0.2]
ad = ad_progression(ad, tf, d, dlimit, thta)

# plotting as a function of t
plt.title('Plot for n=' + str(n) + ' with a=0')
for i in range(n):
        this_label = 'Y' + str(i + 1)
        plt.plot(sol.t, sol.y[i], label=this_label)
plt.xlabel('Time (days)')
plt.ylabel('Mass of Proteins (nanograms)')
plt.legend()

filepath = 'figures/'
os.makedirs(os.path.dirname(filepath), exist_ok=True)
plt.savefig(filepath + 'n=' + str(n))
plt.show()