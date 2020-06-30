import numpy as np
import scipy as sp
from scipy import optimize
import matplotlib.pyplot as plt
from math import *


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
    (T_c, S, T) = data
    return brillouin(x,S) - x * (S+1) * T / (3*S * T_c)


def sol_equation(T):
    T_c = 627
    S = 1
    return sp.optimize.fsolve(equation, 1000, args = (T_c, S, T))


temperature = np.arange(0.01, 627, 1)

liste = []

for T in temperature:
    liste.append(sol_equation(T)[0])


resultat = np.array(liste)
resultat = brillouin(resultat,1)


plt.plot(temperature/627, resultat)
plt.ylabel("M/M_0")
plt.xlabel("T/T_c")
plt.gca()
plt.show()

