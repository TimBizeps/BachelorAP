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

U6, I6 = np.genfromtxt('messwerte4.txt', unpack=True)
I6A = I6/1000000000
U6tat = U6 - 1000000*I6A #1MOhm Innenwiderstand
I6log = np.log(I6A)

e = const.e
k = const.k

def f(x, a, b):
    return a*x+b

params, covariance = curve_fit(f, U6tat, I6log)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

a = ufloat(params[0], errors[0])

T = -e/(k*a)

np.savetxt("Parameter2.txt", np.column_stack([params, errors]))

print('T =', noms(T), '±', stds(T))#in Kelvin

x_plot = np.linspace(-0.05, 1)

plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.plot(U6tat, I6log, 'rx', label='Messwerte', linewidth=1)
plt.xlabel(r'$ln \left( \frac{U}{\si{\volt}} \right)$')
plt.ylabel(r'$ln \left( \frac{I}{\si{\nano\ampere}} \right)$')
plt.xlim(-0.05, 1)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot3.pdf")
