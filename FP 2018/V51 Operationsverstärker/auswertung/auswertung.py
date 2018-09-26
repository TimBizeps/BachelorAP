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

f, Ua = np.genfromtxt('daten1.txt', unpack=True)
f2, Ua2 = np.genfromtxt('daten12.txt', unpack=True)

#für Fit alle Werte f und Vstrich logarithmieren und in nomalen Plot.

Ue = 260

R1 = 468
RN = 469

Vstricheff = Ua/Ue
Vstricheff2 = Ua2/Ue

Vstrich = Ua[0]/Ue

V = 1/((1/Vstrich)-(R1/RN))

print('V = ', V)

Vstrbvg = Vstrich/np.sqrt(2)

print('Vstrich = ', Vstrich)
#print('Vstricheff = ', Vstricheff)
logVstreff = np.log(Vstricheff)
logf = np.log(f)
logVstreff2 = np.log(Vstricheff2)
logf2 = np.log(f2)
logVstrich = np.log(Vstrich)
logVstrbvg = np.log(Vstrbvg)
#print('logVstreff = ', logVstreff)
#print('logf = ', logf)
#np.savetxt("logdaten1.txt", np.column_stack([logf,logVstreff]))
#np.savetxt("logdaten12.txt", np.column_stack([logf2,logVstreff2]))

print('logVstrbvg = ', logVstrbvg)

def linear(x, m, b):
    return m*x+b

def umkehr(y, m, b):
    return (y-b)/m
#fit mit linearer Funktion für die entsprechenden Werte
params, cov = curve_fit(linear, logf2, logVstreff2)

errors = np.sqrt(np.diag(cov))

m = params[0]
m_err = errors[0]

b = params[1]
b_err = errors[1]

mmit = ufloat(m, m_err)
bmit = ufloat(b, b_err)
print('m = ', mmit)
print('b = ', bmit)
#umkehrfunktion mit log Vstrbvg aufrufen => v'g
logvg = umkehr(logVstrbvg, mmit, bmit)
#print('logvg = ', logvg)
vg = unp.exp(logvg)
print('vg = ', vg)

Vstrmalvg = Vstrich*vg

print('Vstrmalvg = ', Vstrmalvg)

l = np.linspace(6, 8, 1000)

plt.plot(logf, logVstreff, 'rx', label='doppelt logarithmische Messwerte')
plt.plot(logf2, logVstreff2, 'kx', label='für den Fit verwendete Messwerte')
plt.plot(l, linear(l,m,b), 'r-', label='Fit')
plt.xlabel(r'$\log \nu$')
plt.ylabel(r'$\log V^\prime_\text{eff}$')
#plt.xlim(,)
#plt.ylim(100, 400)
plt.axhline(y=logVstrbvg, xmin=0, xmax=1, color='b', ls='-', label=r"$\log\frac{V^\prime}{\sqrt{2}}$", linewidth=1)
plt.axhline(y=logVstrich, xmin=0, xmax=1, color='b', ls='--', label=r"$\log V^\prime$", linewidth=1)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot.pdf")

#f_g etwa 1150 kHz
