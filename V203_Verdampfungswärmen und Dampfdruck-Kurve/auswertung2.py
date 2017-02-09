import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

y, z = np.genfromtxt('messwerte2.txt', unpack=True)

z = z+273.15

R = 8.314 #joule/(mol*K)
A = 0.9 #joule*m³/mol²

def f(z, a, b, c, d, e):
    return a * z**4 + b * z**3 + c * z**2 + d * z + e

def g(z, h, i, j, k, l):
    return ((R*z/2)+np.sqrt((R*z/2)**2-A*(h * z**4 + i * z**3 + j * z**2 + k * z + l)))*(4*h * z**4 + 3*i * z**3 + 2*j * z**2 + k * z)/(h * z**4 + i * z**3 + j * z**2 + k * z + l)

def m(z, n, o, p, q, r):
    return ((R*z/2)-np.sqrt((R*z/2)**2-A*(n * z**4 + o * z**3 + p * z**2 + q * z + r)))*(4*n * z**4 + 3*o * z**3 + 2*p * z**2 + q * z)/(n * z**4 + o * z**3 + p * z**2 + q * z + r)

params, covariance = curve_fit(f, z, y)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])
print('d =', params[3], '±', errors[3])
print('e =', params[4], '±', errors[4])

np.savetxt("WertefürL(T).txt", np.column_stack([g(z, *params), m(z, *params)]))

plt.plot(z, m(z, *params), 'rx', label="Fall 2: L(T) bei Subtraktion der Wurzel")
plt.legend(loc="best")
plt.xlabel(r'$T \,/\, K$')
plt.ylabel(r'$L \,/\, Joule/mol$')
plt.savefig("Plot4.pdf")
