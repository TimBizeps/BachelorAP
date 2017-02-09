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

x, y = np.genfromtxt('daten2.txt', unpack=True)

R1 = ufloat(48.1, 0.1)
L = ufloat(0.01011, 0.00003)
C = ufloat(0.000000002098, 0.000000000006)

def f(t, a, b):
    return a*np.exp(-b*t)

params, covariance = curve_fit(f, x, y)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

Exponent = ufloat(params[1], errors[1])

x_plot = np.linspace(0, 0.000245)
plt.plot(x, y, 'rx', label="Messdaten")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Ausgleichskurve', linewidth=1)
plt.legend(loc="best")
plt.xticks([0, 0.00005, 0.00010, 0.00015, 0.00020, 0.00025],
        [0, 50, 100, 150, 200, 250])
plt.xlabel(r'Zeit $t$ in $\si{\micro\second}$')
plt.ylabel(r'Amplitude in $\si{\volt}$')
plt.ylim(15, 90)
plt.tight_layout()
plt.savefig("Plot.pdf")

#print('Exponent = ', Exponent.n, '±', Exponent.s)
Reff = Exponent*2*L
print('Reff = ', Reff.n, '±', Reff.s)
Tex = 1/Exponent
print('Tex = ', Tex.n, '±', Tex.s)#179pm6micros
Rap = 2*unp.sqrt(L/C)
print('Rap = ', Rap.n, '±', Rap.s)
