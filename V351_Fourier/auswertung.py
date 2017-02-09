import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit

n, U = np.genfromtxt('messwerte.txt', unpack=True)
Urel = U/U[0]
logn = np.log(n)
logU = np.log(Urel)

plt.plot(logn, logU, 'rx', label='Messwerte')
plt.xlabel(r'$\mathrm{\log \left( n \right)}$')
#plt.xscale('log')
plt.ylabel(r'$\mathrm{\log \left( \frac{U}{U_0} \right)}$')
#plt.yscale('log')
#plt.ylim(0.02, 10)

def t(x, a, b):
    return a*x+b

params, covariance = curve_fit(t, logn, logU)

errors= np.sqrt(np.diag(covariance))
print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

x_plot = np.linspace(0.615, 12)
# plt.plot(d, U, 'r.', label= 'Messwerte')
plt.plot(np.log(x_plot), t(np.log(x_plot), *params), 'b-', label='Fit')

plt.tight_layout()
plt.savefig('Sägezahn.pdf')
# plt.show()
