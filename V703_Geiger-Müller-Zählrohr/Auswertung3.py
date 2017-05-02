import matplotlib as mpl
mpl.use('pgf')
import numpy as np
from scipy.constants import e
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

x, y, z = np.genfromtxt('Daten3.txt', unpack=True)

zm = np.mean(z)
zstm = np.std(z, ddof=1) / np.sqrt(len(z))

print('Erholungszeit:', zm, '±', zstm)

n11 = 19262
n1 = n11/60
n1f = np.sqrt(n11)
n21 = 513
n2 = n21/60
n2f = np.sqrt(n21)
n121 = 19539
n12 = n121/60
n12f = np.sqrt(n121)

N1 = ufloat(n11, n1f)
N2 = ufloat(n21, n2f)
N12 = ufloat(n121, n12f)

T = 1/N12 + unp.sqrt((1/(N12**2))-((1/(N2*N12))+(1/(N1*N12))-(1/(N1*N2))))
T2 = (N1 + N2 - N12)/(2*N1*N2)

print('T =', T.n, '±', T.s)
print('T2 =', T2.n, '±', T2.s)
