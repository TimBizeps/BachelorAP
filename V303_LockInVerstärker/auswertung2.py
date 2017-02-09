import matplotlib as mpl
mpl.use('pgf')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
mpl.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'lualatex',
    'pgf.preamble': r'\usepackage{unicode-math}\usepackage{siunitx}'
})

x, y = np.genfromtxt('messwerte2.txt', unpack=True)

x=x-4.9
x=x/100

def f(x, a, b):
    return (a/x**b)

params, covariance = curve_fit(f, x, y)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

x_plot=np.linspace(0.05, 1.50)
plt.plot(x, y, 'rx', label="Messwerte")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit-Funktion', linewidth=1)
#plt.yscale('log')
plt.legend(loc="best")
plt.xlabel(r'Distanz zur Photodiode \,/\, $\si{\meter}$')
plt.ylabel(r'Ausgangsspannung \,/\, $\si{\volt}$')
plt.title('Lichtintensitäts Messung')
plt.tight_layout()
plt.savefig("Plot3.pdf")
