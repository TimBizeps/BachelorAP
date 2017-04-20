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

n1 = 19262
n1f = np.sqrt(n1)
n2 = 513
n2f = np.sqrt(n2)
n12 = 19539
n12f = np.sqrt(n12)

N1 = ufloat(n1, n1f)
N2 = ufloat(n2, n2f)
N12 = ufloat(n12, n12f)

T = 1/N12 + unp.sqrt((1/(N12**2))-((1/(N2*N12))+(1/(N1*N12))-(1/(N1*N2))))

print('T =', T.n, '±', T.s)
