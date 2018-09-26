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

x, y = np.genfromtxt('h.txt', unpack=True)
x2, y2 = np.genfromtxt('h2.txt', unpack=True)

def funktion(x,a,b,c,d):
    return a*np.exp((-1*(x-b))/(20*c))*np.sin((x-b)/c)+d

params, cov = curve_fit(funktion, x2, y2, p0=(0.2, 0, 0.000224, -0.02))

errors = np.sqrt(np.diag(cov))

a = params[0]
a_err = errors[0]

b = params[1]
b_err = errors[1]

c = params[2]
c_err = errors[2]

d = params[3]
d_err = errors[3]

print('a = ', a ,'+-', a_err)
print('b = ', b ,'+-', b_err)
print('c = ', c ,'+-', c_err)
print('d = ', d ,'+-', d_err)

l = np.linspace(0, 0.045, 5000)

plt.plot(x, y, 'rx', label='Messwerte, die nicht verwendet werden')
plt.plot(x2, y2, 'kx', label='Messwerte für den Fit')
plt.plot(l, funktion(l,a,b,c,d), 'b', label='Fit')
plt.xlabel(r'Zeit $t$ in $\si{\second}$')
plt.ylabel(r'Ausgangsspannung $U_A$ in $\si{\volt}$')
#plt.xlim(,)
#plt.ylim(100, 500)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot5.pdf")

cmit = ufloat(c, c_err)
ctau = 20*cmit
print('ctau = ', ctau)

k = np.array([215, 236])
meank = np.mean(k)
fk = np.std(k, ddof=1)/np.sqrt(len(k))
kmit = ufloat(meank, fk)
print('kmit = ', kmit)
kmit = kmit * 10**(-10)
konstantec = kmit * 9960
tau = 20*konstantec
print('kmit = ', kmit)
print('konstantec = ', konstantec)
print('tau = ', tau)
