import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit

d, U = np.genfromtxt('messwerte2.txt', unpack=True)

d=d-4.9
d=d/100

logd = np.log(d)
logU = np.log(U)

plt.plot(logd, logU, 'k.', label='$Abstandsverhältnis der Lichtintensität$')
plt.xlabel(r'$\mathrm{log \left( \left| d \right| \right)}$')
#plt.xscale('log')
plt.ylabel(r'$\mathrm{log \left( \left| U \right| \right)}$')
#plt.yscale('log')
#plt.ylim(0.02, 10)

def f(x, b, a):
    return b*x+a

params, cov = curve_fit(f, logd, logU)

errors= np.sqrt(np.diag(cov))
print('			Parameter für Data3 gefittet:	')
print('			b=', params[0], '\pm', errors[0])
print('			a=', params[1], '\pm', errors[1])

x_plot = np.linspace(0.05, 1.50)
# plt.plot(d, U, 'r.', label= 'Messwerte')
plt.plot(np.log(x_plot), f(np.log(x_plot), *params), 'b-', label='Fit')

plt.tight_layout()
plt.savefig('licht-messung.pdf')
# plt.show()
