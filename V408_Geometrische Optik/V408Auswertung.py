import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

g, b_1 = np.loadtxt("V408c.txt", unpack=True) # lädt Messdaten aus txt in Variablen g und b_1

for i in range(0,len(g)):
	plt.plot((g[i], 0), (0, b_1[i]), 'k-')

plt.plot((g[0], 0), (0, b_1[0]), 'k-', label='Verbindungslinie zwischen g und b')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Gegenstandsweite}\; g/\mathrm{cm}$')
plt.ylabel(r'$\mathrm{Bildweite}\; b/\mathrm{cm}$')
plt.tight_layout()

plt.savefig("sam1.pdf")
plt.close()
