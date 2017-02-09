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
m, f = np.genfromtxt('Daten2.txt', unpack=True)

def g(x, a, b):
    return a*np.log(x) + b

def g2(y, a, b):
    return unp.exp((y-b)/a)

params, covariance = curve_fit(g, f, m)
errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

#x_plot = np.linspace(0, 3)
#plt.plot(phase, kreisfreq, 'rx', label="Messdaten")
#plt.plot(x_plot, omega(x_plot, 0.00175, 0.000000022), 'b-', label='Theoriekurve', linewidth=1)
#plt.xscale('log')
#plt.legend(loc="best")
#plt.tight_layout()
#plt.savefig("Plot.pdf")
a = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print(g2(15.9, *params))
print(np.sqrt((2/0.00175)*(0.000000022+0.00000000939)/(0.000000022*0.00000000939))/(2*np.pi))
