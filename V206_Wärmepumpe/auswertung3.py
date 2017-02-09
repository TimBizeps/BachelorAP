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

v, x, n, o = np.genfromtxt('daten2.txt', unpack=True)
n = n + 273.15
o = o + 273.15
p0 = 1.013
x = x/p0
def R(o, h, i):
    return h * o + i

params3, covariance3 = curve_fit(R, 1/o, np.log(x))

errors3 = np.sqrt(np.diag(covariance3))

print('a =', params3[0], '±', errors3[0])
print('b =', params3[1], '±', errors3[1])

x_plot = np.linspace(0.00335, 0.0037)
plt.plot(1/o, np.log(x), 'rx', label="Messwerte")
plt.plot(x_plot, R(x_plot, *params3), 'b-', label='Ausgleichsgerade', linewidth=1)
plt.legend(loc="best")
plt.xlabel(r'$\frac{1}{T} \,\,/\,\, \frac{1}{K}$')
plt.ylabel(r'$ln(\frac{p}{p_0})$')
plt.tight_layout()
plt.savefig("Plot3.pdf")
