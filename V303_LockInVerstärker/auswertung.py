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

x, y, z = np.genfromtxt('messwerte.txt', unpack=True)

x=x*np.pi

def f(x, a, b):
    return a*np.cos(x+b)

params, covariance = curve_fit(f, x, y)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1]+2*np.pi, '±', errors[1])

x_plot=np.linspace(0, 6.9)
plt.plot(x, y, 'rx', label="Messwerte")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit-Funktion', linewidth=1)
plt.legend(loc="best")
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
        [r"$0$", r"$\frac{1}{2} \pi$", r"$\pi$", r"$\frac{3}{2} \pi$", r"$2 \pi$"])
plt.xlabel(r'Phase $\phi$ im Bogenmaß')
plt.ylabel(r'Ausgangsspannung \,/\, $\si{\volt}$')
plt.title('Inaktiver Noise Generator')
plt.tight_layout()
plt.savefig("Plot.pdf")
