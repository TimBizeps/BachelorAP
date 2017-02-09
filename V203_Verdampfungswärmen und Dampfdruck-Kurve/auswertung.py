import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

w, x = np.genfromtxt('messwerte.txt', unpack=True)
p0 = 1.017
x = x+273.15
w = w/p0
def f(x, a, b):
    return a * x + b

params, covariance = curve_fit(f, 1/x, np.log(w))

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

x_plot = np.linspace(0.00261, 0.00345)
plt.plot(1/(x), np.log(w), 'rx', label="Messwerte bis 1 Bar")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Ausgleichsgerade', linewidth=1)
plt.legend(loc="best")
plt.xlabel(r'$\frac{1}{T} \,\,/\,\, \frac{1}{K}$')
plt.ylabel(r'$ln(\frac{p}{p_0})$')
plt.savefig("Plot.pdf")

R = 8.314 #joule/(mol*K)
Steigung = ufloat(params[0], errors[0])
T = 373 #K
Na = 6.022 * 10**23

La = R * T #L_a = W = p*V = R*T
L = -Steigung*R
print('L =', L)
print('La =', La)

Li = L - La
print('Li =', Li.n, '±', Li.s)
L_i = Li/Na
L_i = L_i/(1.602 * 10**-19)
print('Li =', L_i.n, '±', L_i.s)
