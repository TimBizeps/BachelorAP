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

f, p1, p2, h1, h2 = np.genfromtxt('daten.txt', unpack=True)
F, P1, P2, H1, H2 = np.genfromtxt('daten5.txt', unpack=True)
i1 = p1*0.1
i2 = p2*0.1
i12 = h1*0.3
i22 = h2*0.3
I1 = P1*0.1
I2 = P2*0.1
I12 = H1*0.3
I22 = H2*0.3
mu0 = const.mu_0
B1 = mu0*(8/np.sqrt(125))*(11/0.1639)*i1 + mu0*(8/np.sqrt(125))*(154/0.1579)*i12
B2 = mu0*(8/np.sqrt(125))*(11/0.1639)*i2 + mu0*(8/np.sqrt(125))*(154/0.1579)*i22
B12 = mu0*(8/np.sqrt(125))*(11/0.1639)*I1 + mu0*(8/np.sqrt(125))*(154/0.1579)*I12
B22 = mu0*(8/np.sqrt(125))*(11/0.1639)*I2 + mu0*(8/np.sqrt(125))*(154/0.1579)*I22

print('B1 =', B1)
print('B2 =', B2)

def l(x, a, b):
    return a*x+b

params, covariance = curve_fit(l, F, B12)
errors = np.sqrt(np.diag(covariance))
params2, covariance2 = curve_fit(l, F, B22)
errors2 = np.sqrt(np.diag(covariance2))

print('a =', params[0], '±', errors[0]) #Tesla/kHz
print('b =', params[1], '±', errors[1])
print('a2 =', params2[0], '±', errors2[0])
print('b2 =', params2[1], '±', errors2[1])

a1 = ufloat(params[0], errors[0])
a2 = ufloat(params2[0], errors2[0])
h = const.h
muB = 9.274009994*(10**(-24))
print('h =', h)

gf2 = (h/(a2*muB))*1000   #aus h*f = U = gF*muB*B
gf1 = (h/(a1*muB))*1000   #=> B = (h*f)/(gF*muB)
                          #=> a = h/(gF*muB)
                          #=> gF = h/(a*muB)
                          #*1000 wegen Steigung in T/kHz
print('gf1 =', gf1)
print('gf2 =', gf2)

#weiter mit Kapitel 3 der Anleitung => Kernspin
#j = 1/2, s = 1/2, l = 0 wegen Alkali-Atom
#gj = 3.0023*j(j+1)+1.0023(s(s+1)-l(l+1))/2j(j+1)
#gf = gj*(f(f+1)+j(j+1)-i(i+1))/2f(f+1)
#gf = gj*((i+j)*((i+j)+1)+j(j+1)-i(i+1))/(2*(i+j)*((i+j)+1))

gj = (3.0023*((1/2)*((1/2)+1))+1.0023*((1/2)*((1/2)+1)))/(2*((1/2)*((1/2)+1)))
print('gj =', gj)
gftheo1 = gj * (6+(3/4)-(15/4))/12
gftheo2 = gj * (12+(3/4)-(35/4))/24
print('gftheo1 =', gftheo1)
print('gftheo2 =', gftheo2)

I1 = ((gj/gf1)-1)/2
I2 = ((gj/gf2)-1)/2
print('I1 =', I1) #~2  beide etwa 0.5 zu groß 3/2 wäre Lit. 87 Rb
print('I2 =', I2) #~3  beide etwa 0.5 zu groß 5/2 wäre Lit. 85 Rb

#Quadr. ,,, MF = 2 bzw. 3

zee1 = gf1 * muB * l(100, *params)
zee2 = gf2 * muB * l(100, *params)
zee12 = gf1 * muB * l(1000, *params2)
zee22 = gf2 * muB * l(1000, *params2)
quadzee1 = gf1**2 * muB**2 * l(100, *params)**2 * (-3)/(4.53*10**(-24))
quadzee2 = gf2**2 * muB**2 * l(100, *params)**2 * (-5)/(2.01*10**(-24))
quadzee12 = gf1**2 * muB**2 * l(1000, *params2)**2 * (-3)/(4.53*10**(-24))
quadzee22 = gf2**2 * muB**2 * l(1000, *params2)**2 * (-5)/(2.01*10**(-24))
print('Bfeld1 =', l(100, *params))
print('Bfeld2 =', l(1000, *params2))
print('zee1 =', zee1)
print('zee12 =', zee12)
print('zee2 =', zee2)
print('zee22 =', zee22)
print('quadzee1 =', quadzee1)
print('quadzee12 =', quadzee12)
print('quadzee2 =', quadzee2)
print('quadzee22 =', quadzee22)

rb87 = (2.4/7.5)*100 # in %
rb85 = (5.1/7.5)*100 # in %
litrb85 = 72.17
litrb87 = 27.83
abwrb87 = ((rb87/litrb87)-1)*100 # in %
abwrb85 = ((rb85/litrb85)-1)*100 # in %
print('rb87 =', rb87)
print('litrb87 =', litrb87)
print('abwrb87 =', abwrb87)
print('rb85 =', rb85)
print('litrb85 =', litrb85)
print('abwrb85 =', abwrb85)

x_plot = np.linspace(-10,1050)

plt.plot(f, B1*1000000, 'rx', label='Messwerte')
plt.plot(f, B2*1000000, 'bx', label='Messwerte 2')
plt.plot(x_plot, l(x_plot, *params)*1000000, 'g--', label='Fit1', linewidth=1)
plt.plot(x_plot, l(x_plot, *params2)*1000000, 'y--', label='Fit2', linewidth=1)
plt.legend(loc="best")
plt.xlabel(r'Frequenz in $\si{\kilo\hertz}$')
plt.ylabel(r'$B_m + B_{Erde}$ in $\si{\micro\tesla}$')
plt.tight_layout
plt.savefig("Plot1.pdf")
