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
L = 0.023954
C = 0.0000000007932
Csp = 0.000000000028
Ck2 = ufloat(0.000000000997, 0.000000000997*0.03)

Ck, N, F = np.genfromtxt('daten.txt', unpack=True)

Ck = unp.uarray(Ck, 0.03*Ck)
print('Ck =', noms(Ck), '±', stds(Ck))
Ck = Ck*10**(-9)
n = unp.uarray(N, F)
print('Anzahl =', noms(n), '±', stds(n))
nup = (1/(2*np.pi*np.sqrt(L*(C+Csp))))/1000
print('nu+ =', nup)
num = (1/(2*np.pi*unp.sqrt(L*((1/((1/C)+(2/Ck)))+Csp))))/1000
print('nu-theo =', noms(num), '±', stds(num))
nuosz = (1/2)*(nup+num)
nusch = num-nup
nt = nuosz/nusch
print('nt =', noms(nt), '±', stds(nt))

nup = np.array([nup, nup, nup, nup, nup, nup, nup])

m = 66-25 #kHz/s
b = 25 #kHz

def f(t, m, b):
    return m * t + b

t11 = ufloat(0.252, 0.004)
t12 = ufloat(0.740, 0.004)
t21 = ufloat(0.260, 0.004)
t22 = ufloat(0.512, 0.004)
t31 = ufloat(0.256, 0.004)
t32 = ufloat(0.460, 0.004)
t41 = ufloat(0.256, 0.004)
t42 = ufloat(0.384, 0.004)
t51 = ufloat(0.256, 0.004)
t52 = ufloat(0.344, 0.004)
t61 = ufloat(0.260, 0.004)
t62 = ufloat(0.336, 0.004)
t71 = ufloat(0.256, 0.004)
t72 = ufloat(0.316, 0.004)
t81 = ufloat(0.256, 0.004)
t82 = ufloat(0.304, 0.004)

nup1 = f(t11,m,b)#0.997nF
num1 = f(t12,m,b)
nup2 = f(t21,m,b)#2.19
num2 = f(t22,m,b)
nup3 = f(t31,m,b)#2.86
num3 = f(t32,m,b)
nup4 = f(t41,m,b)#4.74
num4 = f(t42,m,b)
nup5 = f(t51,m,b)#6.86
num5 = f(t52,m,b)
nup6 = f(t61,m,b)#8.18
num6 = f(t62,m,b)
nup7 = f(t71,m,b)#9.99
num7 = f(t72,m,b)
nup8 = f(t81,m,b)#12
num8 = f(t82,m,b)

numtheo2 = 1/(2*np.pi*unp.sqrt(L*((1/((1/C)+(2/Ck2)))+Csp)))
print('nu-theo2 =', numtheo2.n, '±', numtheo2.s)

print('nu+1 =', nup1.n, '±', nup1.s)
print('nu-1 =', num1.n, '±', num1.s)
print('nu+2 =', nup2.n, '±', nup2.s)
print('nu-2 =', num2.n, '±', num2.s)
print('nu+3 =', nup3.n, '±', nup3.s)
print('nu-3 =', num3.n, '±', num3.s)
print('nu+4 =', nup4.n, '±', nup4.s)
print('nu-4 =', num4.n, '±', num4.s)
print('nu+5 =', nup5.n, '±', nup5.s)
print('nu-5 =', num5.n, '±', num5.s)
print('nu+6 =', nup6.n, '±', nup6.s)
print('nu-6 =', num6.n, '±', num6.s)
print('nu+7 =', nup7.n, '±', nup7.s)
print('nu-7 =', num7.n, '±', num7.s)
print('nu+8 =', nup8.n, '±', nup8.s)
print('nu-8 =', num8.n, '±', num8.s)

nupg2n = np.array([nup8.n, nup7.n, nup6.n, nup5.n, nup4.n, nup3.n, nup2.n])
numg2n = np.array([num8.n, num7.n, num6.n, num5.n, num4.n, num3.n, num2.n])
nupg2s = np.array([nup8.s, nup7.s, nup6.s, nup5.s, nup4.s, nup3.s, nup2.s])
numg2s = np.array([num8.s, num7.s, num6.s, num5.s, num4.s, num3.s, num2.s])

kc, nupg, numg, er = np.genfromtxt('daten2.txt', unpack=True)

plt.errorbar(noms(Ck), nup, xerr=stds(Ck), yerr=0, fmt='bo', label=r"$\nu_+$ theoretisch")
plt.errorbar(noms(Ck), noms(num), xerr=stds(Ck), yerr=stds(num), fmt='ko', label=r"$\nu_-$ theoretisch")
plt.errorbar(noms(Ck), nupg, xerr=stds(Ck), yerr=er, fmt='co', label=r"$\nu_+$ direkt gemessen")
plt.errorbar(noms(Ck), numg, xerr=stds(Ck), yerr=er, fmt='ro', label=r"$\nu_-$ direkt gemessen")
plt.errorbar(noms(Ck), nupg2n, xerr=stds(Ck), yerr=nupg2s, fmt='go', label=r"$\nu_+$ per Sweep gemessen")
plt.errorbar(noms(Ck), numg2n, xerr=stds(Ck), yerr=numg2s, fmt='yo', label=r"$\nu_-$ per Sweep gemessen")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot.pdf")
