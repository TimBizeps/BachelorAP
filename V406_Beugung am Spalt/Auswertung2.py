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

y, x = np.genfromtxt('Messwerte2.txt', unpack=True)

ys = y-0.00125

def I(x, a, b, c, d):
    return a * (np.cos((np.pi*d*np.sin((x-c)/1240))/0.000635))**2 * (0.000635/(np.pi*b*np.sin((x-c)/1240)))**2 * (np.sin((np.pi*b*np.sin((x-c)/1240))/0.000635))**2

params, covariance = curve_fit(I, x, ys, p0=[1.8, 0.1, 27, 0.4])

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])
print('d =', params[3], '±', errors[3])

x_plot = np.linspace(15, 40.5, 1000)

np.savetxt("Parameter2.txt", np.column_stack([params, errors]))

plt.plot(x_plot, I(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.plot(x, ys, 'rx', label='Messwerte', linewidth=1)
plt.xlabel(r'$x$ / $\si{\milli\meter}$')
plt.ylabel(r'$I$ / $\si{\micro\ampere}$')
plt.xlim(15, 40.5)
plt.ylim(0, 2.4)
plt.grid()
plt.legend(loc="best")
plt.title('Beugung am Doppelspalt 1')
plt.tight_layout()
plt.savefig("Plot2.pdf")
