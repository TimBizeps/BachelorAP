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

B, I = np.genfromtxt('daten.txt', unpack=True)

def linear(x,a,b):
    return a*x+b

params, cov = curve_fit(linear, I, B)

m = params[0]
m_err = np.sqrt(cov[0][0])

b = params[1]
b_err = np.sqrt(cov[1][1])

c = np.linspace(0, 20)

plt.plot(I,B, 'rx', label='Messwerte')
plt.plot(c, linear(c,m,b), label='Regression')
plt.xlabel(r'Stromstärke $I$ in $\si{\ampere}$')
plt.ylabel(r'Magnetfeld $B$ in $\si{\milli\tesla}$')
plt.xlim(0, 20)
plt.ylim(0, 1050)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot.pdf")

print('b ', b)
print('b_err ', b_err)
print('m ', m)
print('m_err ', m_err)

bmit = ufloat(b, b_err)
mmit = ufloat(m, m_err)
Ir = 11.25
I1 = 9.75
I2 = 13.4
I3 = 15.25
I4 = 16
I5 = 19.4

Br = linear(Ir, mmit, bmit)/1000 #Tesla
B1 = linear(I1, mmit, bmit)/1000 #Tesla
B2 = linear(I2, mmit, bmit)/1000 #Tesla
B3 = linear(I3, mmit, bmit)/1000 #Tesla
B4 = linear(I4, mmit, bmit)/1000 #Tesla
B5 = linear(I5, mmit, bmit)/1000 #Tesla

print('Br ', Br)
print('B1 ', B1)
print('B2 ', B2)
print('B3 ', B3)
print('B4 ', B4)
print('B5 ', B5)

DeltaSr, deltas = np.genfromtxt('daten3.txt', unpack=True)

quoti = deltas/DeltaSr

DeltaLambDr = 4.89*10**(-11)
DeltaLambDb = 2.7*10**(-11)

print('rot ', quoti)

DeltaSb, deltas1, deltas2, deltas31, deltas32, deltas4 = np.genfromtxt('daten4.txt', unpack=True)

quoti1 = deltas1/DeltaSb
quoti2 = deltas2/DeltaSb
quoti31 = deltas31/DeltaSb
quoti32 = deltas32/DeltaSb
quoti4 = deltas4/DeltaSb

print('b1 ', quoti1)
print('b2 ', quoti2)
print('b31 ', quoti31)
print('b32 ', quoti32)
print('b4 ', quoti4)

DeltaSb2, deltas5 = np.genfromtxt('daten5.txt', unpack=True)

quoti5 = deltas5/DeltaSb2

print('b5 ', quoti5)

deltaLambr = (1/2)*quoti*DeltaLambDr       #1/2 weil die Aufspaltung nach links
deltaLambb1 = (1/2)*quoti1*DeltaLambDb     #und rechts geht und deltaLamb dann zur mitte gemessenwird
deltaLambb2 = (1/2)*quoti2*DeltaLambDb
deltaLambb31 = (1/2)*quoti31*DeltaLambDb
deltaLambb32 = (1/2)*quoti32*DeltaLambDb
deltaLambb4 = (1/2)*quoti4*DeltaLambDb
deltaLambb5 = (1/2)*quoti5*DeltaLambDb
print('deltaLambr ', deltaLambr)
print('deltaLambb1 ', deltaLambb1)
print('deltaLambb2 ', deltaLambb2)
print('deltaLambb31 ', deltaLambb31)
print('deltaLambb32 ', deltaLambb32)
print('deltaLambb4 ', deltaLambb4)
print('deltaLambb5 ', deltaLambb5)
deltaLambb3 = np.concatenate((deltaLambb31, deltaLambb32), axis=None)
print('deltaLambb3 ', deltaLambb3)
meandeltalar = np.mean(deltaLambr)
fdeltalar = np.std(deltaLambr, ddof=1)/np.sqrt(len(deltaLambr))
meandeltalab1 = np.mean(deltaLambb1)
fdeltalab1 = np.std(deltaLambb1, ddof=1)/np.sqrt(len(deltaLambb1))
meandeltalab2 = np.mean(deltaLambb2)
fdeltalab2 = np.std(deltaLambb2, ddof=1)/np.sqrt(len(deltaLambb2))
meandeltalab3 = np.mean(deltaLambb3)
fdeltalab3 = np.std(deltaLambb3, ddof=1)/np.sqrt(len(deltaLambb3))
meandeltalab4 = np.mean(deltaLambb4)
fdeltalab4 = np.std(deltaLambb4, ddof=1)/np.sqrt(len(deltaLambb4))
meandeltalab5 = np.mean(deltaLambb5)
fdeltalab5 = np.std(deltaLambb5, ddof=1)/np.sqrt(len(deltaLambb5))

deltaLaR = ufloat(meandeltalar, fdeltalar)
deltaLaB1 = ufloat(meandeltalab1, fdeltalab1)
deltaLaB2 = ufloat(meandeltalab2, fdeltalab2)
deltaLaB3 = ufloat(meandeltalab3, fdeltalab3)
deltaLaB4 = ufloat(meandeltalab4, fdeltalab4)
deltaLaB5 = ufloat(meandeltalab5, fdeltalab5)

print('deltaLaR ', deltaLaR)
print('deltaLaB1 ', deltaLaB1)
print('deltaLaB2 ', deltaLaB2)
print('deltaLaB3 ', deltaLaB3)
print('deltaLaB4 ', deltaLaB4)
print('deltaLaB5 ', deltaLaB5)

h = const.h
c = const.c
müB = 9.274009994*10**(-24) #J/T

#E = (h*c)/Lambda
#E = gj*müB*B
#=> gj = (h*c)/(Lambda*müB*B)

deltaEr = ((h*c)/(643.8*10**(-9))-(h*c)/(deltaLaR+(643.8*10**(-9))))
deltaEb1 = ((h*c)/(480*10**(-9))-(h*c)/(deltaLaB1+(480*10**(-9))))
deltaEb2 = ((h*c)/(480*10**(-9))-(h*c)/(deltaLaB2+(480*10**(-9))))
deltaEb3 = ((h*c)/(480*10**(-9))-(h*c)/(deltaLaB3+(480*10**(-9))))
deltaEb4 = ((h*c)/(480*10**(-9))-(h*c)/(deltaLaB4+(480*10**(-9))))
deltaEb5 = ((h*c)/(480*10**(-9))-(h*c)/(deltaLaB5+(480*10**(-9))))

EBr = müB*Br
EB1 = müB*B1
EB2 = müB*B2
EB3 = müB*B3
EB4 = müB*B4
EB5 = müB*B5

print('deltaE ', deltaEr)
print(deltaEb1)
print(deltaEb2)
print(deltaEb3)
print(deltaEb4)
print(deltaEb5)

print('EB ', EBr)
print(EB1)
print(EB2)
print(EB3)
print(EB4)
print(EB5)

gjR = deltaEr/EBr
gjB1 = deltaEb1/EB1
gjB2 = deltaEb2/EB2
gjB3 = deltaEb3/EB3
gjB4 = deltaEb4/EB4
gjB5 = deltaEb5/EB5

print('gjR ', gjR)
print('gjB1 ', gjB1)
print('gjB2 ', gjB2)
print('gjB3 ', gjB3)
print('gjB4 ', gjB4)
print('gjB5 ', gjB5)
