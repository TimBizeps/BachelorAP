import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

a, b = np.genfromtxt('DatenRaumT.txt', unpack=True)

t_gr = np.mean(a)
delta_t_gr = 1/np.sqrt(len(a)) * np.std(a, ddof=1)
t_kl = np.mean(b)
delta_t_kl = 1/np.sqrt(len(b)) * np.std(b, ddof=1)

c, d, e, f, g, h, i, j, k, l = np.genfromtxt('DatenTemperaturen.txt', unpack=True)

t_gr_25 = np.mean(c)
delta_t_25 = 1/np.sqrt(len(c)) * np.std(c, ddof=1)
t_gr_30 = np.mean(d)
delta_t_30 = 1/np.sqrt(len(d)) * np.std(d, ddof=1)
t_gr_35 = np.mean(e)
delta_t_35 = 1/np.sqrt(len(e)) * np.std(e, ddof=1)
t_gr_40 = np.mean(f)
delta_t_40 = 1/np.sqrt(len(f)) * np.std(f, ddof=1)
t_gr_45 = np.mean(g)
delta_t_45 = 1/np.sqrt(len(g)) * np.std(g, ddof=1)
t_gr_50 = np.mean(h)
delta_t_50 = 1/np.sqrt(len(h)) * np.std(h, ddof=1)
t_gr_55 = np.mean(i)
delta_t_55 = 1/np.sqrt(len(i)) * np.std(i, ddof=1)
t_gr_60 = np.mean(j)
delta_t_60 = 1/np.sqrt(len(j)) * np.std(j, ddof=1)
t_gr_65 = np.mean(k)
delta_t_65 = 1/np.sqrt(len(k)) * np.std(k, ddof=1)
t_gr_70 = np.mean(l)
delta_t_70 = 1/np.sqrt(len(l)) * np.std(l, ddof=1)

Zeiten = [t_kl, t_gr, t_gr_25, t_gr_30, t_gr_35, t_gr_40, t_gr_45, t_gr_50, t_gr_55, t_gr_60, t_gr_65, t_gr_70]
Fehler = [delta_t_kl, delta_t_gr, delta_t_25, delta_t_30, delta_t_35, delta_t_40, delta_t_45, delta_t_50, delta_t_55, delta_t_60, delta_t_65, delta_t_70]

np.savetxt("DatenPlot.txt", np.column_stack([Zeiten, Fehler]))
