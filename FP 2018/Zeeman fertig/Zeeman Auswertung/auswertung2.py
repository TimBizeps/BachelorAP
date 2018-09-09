import matplotlib as mpl
mpl.use('pgf')
import numpy as np
import scipy.constants as const
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

n, l = np.genfromtxt('daten2.txt', unpack=True)

def f(x,a,b,c):
    return a*np.exp(b*x)+c

params, cov = curve_fit(f, n, l)

a = params[0]
a_err = np.sqrt(cov[0][0])

b = params[1]
b_err = np.sqrt(cov[1][1])

c = params[2]
c_err = np.sqrt(cov[2][2])

werte = np.linspace(0, 20)

plt.plot(n,l, 'rx', label='Messwerte')
plt.plot(werte, f(werte,a,b,c), label='Regression')
plt.xlabel(r'Nummer der Linie')
plt.ylabel(r'Abstand in $\si{\centi\meter}$')
plt.xlim(0, 10)
plt.ylim(0, 2)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot2.pdf")

print('a ', a)
print('a_err ', a_err)
print('b ', b)
print('b_err ', b_err)
print('c ', c)
print('c_err ', c_err)

amit = ufloat(a, a_err)
bmit = ufloat(b, b_err)
cmit = ufloat(c, c_err)
