import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

values = np.array([2.2,
3.2,
1.6,
1.4,
1.3,
1.2,
1.1,
1.0,
1.0,
0.9])
mittelwert = np.mean(values)
varianz = np.var(values)

print(mittelwert)
print(varianz)
