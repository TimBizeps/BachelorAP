import matplotlib as mpl
mpl.use('pgf')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
mpl.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'lualatex',
    'pgf.preamble': r'\usepackage{unicode-math}\usepackage{siunitx}'
})

y, x = np.genfromtxt('Messwerte.txt', unpack=True)
l, k = np.genfromtxt('Messwerte2.txt', unpack=True)

ys = y-0.00115
ls = l-0.00125

def I(x, a, b, c):
    return a * b**2 * (0.000635/(np.pi*b*np.sin((x-c)/1240)))**2 * (np.sin((np.pi*b*np.sin((x-c)/1240))/0.000635))**2

params, covariance = curve_fit(I, x, ys, p0=[0.99, 0.075, 26])

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])

def I2(x, e, f, g, h):
    return e * (np.cos((np.pi*h*np.sin((x-g)/1240))/0.000635))**2 * (0.000635/(np.pi*f*np.sin((x-g)/1240)))**2 * (np.sin((np.pi*f*np.sin((x-g)/1240))/0.000635))**2

Params, Covariance = curve_fit(I2, k, ls, p0=[1.8, 0.1, 27, 0.4])

Errors = np.sqrt(np.diag(Covariance))

print('A =', Params[0], '±', Errors[0])
print('B =', Params[1], '±', Errors[1])
print('C =', Params[2], '±', Errors[2])
print('D =', Params[3], '±', Errors[3])

x_plot = np.linspace(0, 50, 1000)

plt.plot(x_plot, 1.35*I(x_plot, 170.068, 0.0999, 27.0767), 'r-', label='Theoretisch: Einzelspalt', linewidth=1)
plt.plot(x_plot, I2(x_plot, *Params), 'b-', label='Fit des ersten Doppelspalts', linewidth=1)
plt.plot(k, ls, 'gx', label='Messwerte Doppelspalt', linewidth=1)
plt.xlabel(r'$x$ / $\si{\milli\meter}$')
plt.ylabel(r'$I$ / $\si{\micro\ampere}$')
#plt.yscale('log') (wie das Auge es sieht)
plt.xlim(0, 50)
plt.ylim(0, 2.4)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot4.pdf")
