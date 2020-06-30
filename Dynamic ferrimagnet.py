import numpy as np
import scipy as sp
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

tho_0 = 100 * 10^(-15)
C_A = 1/(tho_0)
C_B = 10/(tho_0)
tho_e = 10 * tho_0
T_emax = 500

def rateup(nu, updown, m) :
    if nu == 'A' and updown == 'up' :
        resultat = C_A * (S_A * (S_A + 1) - m * (m + 1)) * 1 / (k_B * T_e[t]) * (
                    n_c * alpha_A * (p * khi_AA * x[0] + (1 - p) * khi_AB * x[1])) / (1 - e ** (
                    1 / (k_B * T_e[t]) * -(n_c * alpha_A * (p * khi_AA * x[0] + (1 - p) * khi_AB * x[1]))))
        return resultat
    if nu == 'A' and updown == 'down' :
        resultat = C_A * (S_A * (S_A + 1) - m * (m - 1)) * 1 / (k_B * T_e[t]) * -(
                    n_c * alpha_A * (p * khi_AA * x[0] + (1 - p) * khi_AB * x[1])) / (1 - e ** (
                    1 / (k_B * T_e[t]) * (n_c * alpha_A * (p * khi_AA * x[0] + (1 - p) * khi_AB * x[1]))))
        return resultat
    if nu == 'B' and updown == 'up':
        resultat = C_B * (S_B * (S_B + 1) - m * (m + 1)) * 1 / (k_B * T_e[t]) * (
                n_c * alpha_B * (p * khi_BA * x[0] + (1 - p) * khi_BB * x[1])) / (1 - e ** (
                1 / (k_B * T_e[t]) * -(n_c * alpha_B * (p * khi_BA * x[0] + (1 - p) * khi_BB * x[1]))))
        return resultat
    if nu == 'B' and updown == 'down':
        resultat = C_B * (S_B * (S_B + 1) - m * (m - 1)) * 1 / (k_B * T_e[t]) * -(
                n_c * alpha_B * (p * khi_BA * x[0] + (1 - p) * khi_BB * x[1])) / (1 - e ** (
                1 / (k_B * T_e[t]) * (n_c * alpha_B * (p * khi_BA * x[0] + (1 - p) * khi_BB * x[1]))))
        return resultat

def T_e(T_0, T_emax,temps):
    T_e = {}
    for t in temps:
        T_e[t] = T_0 + T_emax (1 - e**(-t/(tho_0))) * e**(-t/tho_e) )
    return T_e

