import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

x, y = np.genfromtxt('literaturwerteeta.txt', unpack=True)

def f(x, A, B):
    return A+B*x

params, covariance = curve_fit(f, 1/x, np.log(y))

errors = np.sqrt(np.diag(covariance))

lna = ufloat(params[0], errors[0])
a = unp.exp(lna)

print('A =', params[0], '±', errors[0])
print('B =', params[1], '±', errors[1])

print('a =', unp.nominal_values(a), '±', unp.std_devs(a))

np.savetxt("Parameter.txt", np.column_stack([params, errors]))

x_plot = np.linspace(0.0029, 0.0035)

plt.plot(1/x, np.log(y), 'ro', label='Ergebnisse')
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.legend(loc="best")
plt.xlabel(r'$\frac{1}{T} \,/\, K$')
plt.ylabel(r'$Viskosität \,/\, ln(\eta$)')

plt.savefig("Plot2.pdf")
