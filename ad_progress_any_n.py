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
        deg = 0
        for s in range(0, len(odeys) - 1):
                deg += gammas[s]*odeys[s][len(odeys[s])-1]
        return deg


def ad_progression(ad, endtime, deg, deglimit, theta):
        # updated twice a year
        t = 0.5
        while t <= endtime:
                deg_val = int(t*400)

                t += 0.5
        return ad


n = int(input('input a value for n: '))

# the final time and time step for evaluating
years = 20
tslice = 0.5
curr_slice = 0.0
dt = 0.1

# defining constants
k = 10 ** (-4)
kstar = 5 * 10 ** (-6)
M = 10 ** (-2)
beta = 15

# set gamma influences to be uniform among oligomers
gmma = [0]
for i in range(n-2):
        gmma.append(1/(n-2))
gmma.append(0)
print("gammas: " + str(gmma))

# set ad progression constants
dlimit = 22
thta = 10**(-3)

# initializing a(t) and odes
ad = [0.02]
prev_ad = 0.02
initial = tuple([0] * n)

# initializing solution
tsol = []
ysol = {}

while curr_slice < years:
        # lmda set with current a(t) value
        lmda = 2 * (1 - prev_ad) * (1 + beta * prev_ad)
        # reassigning constants with new lmda
        constants = (k, kstar, lmda, M)
        print("current slice: " + str(curr_slice))
        if curr_slice != 0:
                start_t = curr_slice
                end_slice = curr_slice + tslice + dt
        else:
                start_t = curr_slice
                end_slice = curr_slice + tslice + dt
        print("(" + str(start_t) + "," + str(end_slice) + ")")
        next_sol = integrate.solve_ivp(system, (start_t, end_slice), initial,
                                       t_eval=np.arange(start_t, end_slice, dt), args=constants)
        if len(tsol) != 0:
                print(tsol)
                print(next_sol.t)
                tsol = np.append(tsol, next_sol.t[1:len(next_sol.t)-1])
                for i in range(len(ysol)):
                        print(ysol[i][len(ysol[i])-1])
                        print(next_sol.y[i])
                        ysol[i] = np.append(ysol[i], next_sol.y[i][1:len(next_sol.y[i])-1])
        else:
                tsol = next_sol.t
                for i in range(len(next_sol.y)):
                        ysol[i] = next_sol.y[i]
        # update initial for next slice
        new_init = []
        for i in range(len(ysol)):
                new_init.append(ysol[i][len(ysol[i])-1])
        initial = tuple(new_init)
        curr_slice += tslice
        # calculate degradation
        d = degradation(gmma, ysol)
        # update ad with new_ad
        print("d: " + str(d))
        for i in range(int(tslice/dt)-1):
                ad.append(None)
        if d > dlimit:
                new_ad = prev_ad + theta * (deg[deg_val] - dlimit)
                prev_ad = new_ad
                ad.append(new_ad)
        else:
                ad.append(prev_ad)

ad.pop()
print(tsol)
print(ysol[0])
print("times: " + str(len(tsol)))
print("ad's: " + str(len(ad)))
# plotting ys over time
plt.title('Plot for n=' + str(n) + ' with a=0')
for i in range(n):
        this_label = 'Y' + str(i + 1)
        plt.plot(tsol, ysol[i], label=this_label)
plt.xlabel('Time (days)')
plt.ylabel('Mass of Proteins (nanograms)')
plt.legend()

# plotting degredation values over time
plt.title('AD Progression over Time')
this_label = 'a(t)'
plt.scatter(tsol, ad, label=this_label)
plt.xlabel('Time (years)')
plt.ylabel('Degradation')
plt.legend()

filepath = 'figures/'
os.makedirs(os.path.dirname(filepath), exist_ok=True)
plt.savefig(filepath + 'ad_progress_immobile_n-' + str(n))
plt.show()