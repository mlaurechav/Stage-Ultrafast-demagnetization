import numpy as np
import scipy as sp
from scipy import optimize
import matplotlib.pyplot as plt
from math import *

khi_AA = 1
khi_BB = 1
khi_AB = -1
khi_BA = -0.1

n_c = 10**(29)

alpha_A = 0.01 * 1.6 * 10**(-19) / n_c
alpha_B = 0.1 * 1.6 * 10**(-19) / n_c

k_B = 1.38 * 10**(-23)

S_A = 3
S_B = 1

p = 0.35


def brillouin(x,J) :
    ind = -1
    for j in range(0,np.size(x)-1) :
       if x[j] == 0 :
            ind = j
            x = np.delete(x,ind)
    resultat = (2*J+1)/(2*J)*(1/np.tanh((2*J+1)/(2*J)*x)) - 1/(2*J)*(1/np.tanh(x/(2*J)))
    if ind != -1 :
        resultat = np.insert(resultat, ind, 0)
        x = np.insert(x, ind, 0)
    return resultat

def equation(x, *data):
    (S_A, S_B, k_B, alpha_A, alpha_B, n_c, khi_AA, khi_AB, khi_BA, khi_BB, p, T) = data
    return [x[0] - S_A * brillouin( (S_A * n_c* alpha_A) / (k_B * T) * (p * khi_AA * x[0] + (1-p) * khi_AB * x[1]),S_A) ,
            x[1] - S_B * brillouin( (S_B * n_c* alpha_B)/ (k_B * T) * (p * khi_BA * x[0] + (1-p) * khi_BB * x[1]),S_B)]

def sol_equation(T):
    solution = []
    for T in temperature :
        data = (S_A, S_B, k_B, alpha_A, alpha_B, n_c, khi_AA, khi_AB, khi_BA, khi_BB, p, T)
        solution.append(sp.optimize.fsolve(equation, [100, 100], args = data))
    return solution

temperature = np.arange(0.01, 550, 1)

resultat = sol_equation(temperature)

P_A = []
P_B = []
P_tot = []

for i in range(len(resultat)):
    P_A.append(resultat[i][0]*p)
    P_B.append(resultat[i][1]*(1-p))
    P_tot.append(P_A[i] + P_B[i])


plt.plot(temperature, P_A, P_B)
plt.plot(temperature, P_tot)
plt.xlabel("Temperature")
plt.show()





