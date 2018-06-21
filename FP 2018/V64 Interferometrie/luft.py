import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy
import uncertainties.unumpy as unp
from scipy.stats import stats
import math as math

N = np.genfromtxt('luft.txt', unpack=True)

N = unumpy.uarray(N, 2)

def n(N, lam, L):
    return N * lam / L + 1


lam = 633 * 10**(-9)
L = 0.1

n_array = n(N, lam, L)

n_best = np.mean(n_array)
print(n_array)
print(n_best)

np.savetxt('luft1.txt', np.column_stack([
        unp.nominal_values(N),
        unp.std_devs(N),
        unp.nominal_values(n_array),
        unp.std_devs(n_array),
        ]), delimiter=' & ', newline=r' \\'+'\n',
        fmt='%.0f & %.0f & %.7f & %.7f')
