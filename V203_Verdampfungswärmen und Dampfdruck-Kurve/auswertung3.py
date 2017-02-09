import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

y, z = np.genfromtxt('messwerte2.txt', unpack=True)

z = z+273.15

R = 8.314 #joule/(mol*K)
A = 0.9 #joule*m³/mol²

def f(z, a, b, c, d):
    return a * z**3 + b * z**2 + c * z + d

def g(z, h, i, j, k):
    return ((R*z/2)+np.sqrt((R*z/2)**2-A*(h * z**3 + i * z**2 + j * z + k)))*(3*h * z**3 + 2*i * z**2 + j * z)/(h * z**3 + i * z**2 + j * z + k)

def t(z, n, o, p, q):
    return ((R*z/2)-np.sqrt((R*z/2)**2-A*(n * z**3 + o * z**2 + p * z + q)))*(3*n * z**3 + 2*o * z**2 + p * z)/(n * z**3 + o * z**2 + p * z + q)

params, covariance = curve_fit(f, z, y)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])
print('d =', params[3], '±', errors[3])

plt.plot(z, g(z, *params), 'rx', label="Fall 1: L(T) bei Addition der Wurzel bei Polynom 3. Grades")
plt.legend(loc="best")
plt.xlabel(r'$T \,/\, K$')
plt.ylabel(r'$L \,\,/\,\, Joule/mol$')
plt.savefig("Plot5.pdf")
