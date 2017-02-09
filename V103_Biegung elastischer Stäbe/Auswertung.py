import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

from scipy.optimize import curve_fit

def auswertung(material, querschnitt, einspannung, x, D, d, L, M):

    if einspannung == "einseitig":

        u = L*x**2 - x**3/3
        g = const.g
        F = M*g
#        d = np.mean(k)
#        Δd = np.sqrt(1/(len(k)*(len(k)-1))*sum((d-k)**2))

        if querschnitt == "kreisfoermig":
            I = np.pi/64*d**4
#            ΔI = np.pi/16*d**3*Δd
        if querschnitt == "quadratisch":
            I = d**4/12
#            ΔI = 1/3*d**3*Δd

        def f(x, m, b):
            return m*x + b

        params, cov = curve_fit(f, u, D)
        m  = params[0]
        b  = params[1]
        Δm = np.sqrt(cov[0][0])
        Δb = np.sqrt(cov[1][1])

        E  = F/(2*I*m)
#        ΔE = np.sqrt((F/(2*I**2*m)*ΔI)**2+(F/(2*I*m**2)*Δm)**2)
        ΔE = np.sqrt((F/(2*I*m**2)*Δm)**2)

        t = np.linspace(u.min(), u.max(), 1000)
        plt.plot(u, 1000*D, 'rx', label='Messwerte')
        plt.plot(t, 1000*f(t, m, b), 'k-', label='Regressionsgerade')
        plt.xlim(u.min(), u.max())
        plt.xlabel(r"$(Lx^2 - \frac{x^3}{3})/\mathrm{m}^3$")
        plt.ylabel(r"$D/\mathrm{mm}$")
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig("build/plot_{}_{}_{}.pdf".format(material, querschnitt, einspannung))
        plt.close()

        print(
        """
        ------------------------------------------------------------------------
        Material:                   {}
        Querschnitt:                {}
        Einspannung:                {}
        Durchmesser d:              {} ± {:.5f} mm
        Länge L:                    {} cm
        Masse M:                    {} kg
        Flächenträgheitsmoment I:   {} ± {} mm^4
        Elastizitätsmodul E:   {} ± {} N/m^2
        Steigung m:            {} ± {}
        Ordinatenabschnitt b:  {} ± {}
        ------------------------------------------------------------------------
        """.format(material, querschnitt, einspannung, d*1e3, 0, L*1e2, M, I*1e12, 0, E*1e0, ΔE*1e0, m, Δm, b, Δb))

    if einspannung == "beidseitig":

        x1, x2 = np.array_split(x, 2)
        D1, D2 = np.array_split(D, 2)

        u1 = 3*L**2*x1 - 4*x1**3
        u2 = 4*x2**3 - 12*L*x2**2 + 9*L**2*x2 - L**3
        g  = const.g
        F  = M*g
#        d = np.mean(k)
#        Δd = np.sqrt(1/(len(k)*(len(k)-1))*sum((d-k)**2))

        if querschnitt == "kreisfoermig":
            I = np.pi/64*d**4
#            ΔI = np.pi/16*d**3*Δd
        if querschnitt == "quadratisch":
            I = d**4/12
#            ΔI = 1/3*d**3*Δd

        def f(x, m, b):
            return m*x + b

        params1, cov1 = curve_fit(f, u1, D1)
        params2, cov2 = curve_fit(f, u2, D2)
        m1  = params1[0]
        m2  = params2[0]
        b1  = params1[1]
        b2  = params2[1]
        Δm1 = np.sqrt(cov1[0][0])
        Δm2 = np.sqrt(cov2[0][0])
        Δb1 = np.sqrt(cov1[1][1])
        Δb2 = np.sqrt(cov2[1][1])

        E1  = F/(48*I*m1)
        E2  = F/(48*I*m2)
#        ΔE1 = np.sqrt((F/(48*I**2*m1)*ΔI)**2+(F/(48*I*m1**2)*Δm1)**2)
        ΔE1 = np.sqrt((F/(48*I*m1**2)*Δm1)**2)
#        ΔE2 = np.sqrt((F/(48*I**2*m2)*ΔI)**2+(F/(48*I*m2**2)*Δm2)**2)
        ΔE2 = np.sqrt((F/(48*I*m2**2)*Δm2)**2)

        E  = (E1+E2)/2
        ΔE = np.sqrt(ΔE1**2+ΔE2**2)/2

        t = np.linspace(u1.min(), u1.max(), 1000)
        plt.plot(u1, 1000*D1, 'rx', label='Messwerte')
        plt.plot(t, 1000*f(t, m1, b1), 'k-', label='Regressionsgerade')
        plt.xlim(u1.min(), u1.max())
        plt.xlabel(r"$(3L^2x - 4x^3)/\mathrm{m}^3$")
        plt.ylabel(r"$D/\mathrm{mm}$")
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig("build/plot_{}_{}_{}_1.pdf".format(material, querschnitt, einspannung))
        plt.close()

        t = np.linspace(u2.min(), u2.max(), 1000)
        plt.plot(u2, 1000*D2, 'rx', label='Messwerte')
        plt.plot(t, 1000*f(t, m2, b2), 'k-', label='Regressionsgerade')
        plt.xlim(u2.min(), u2.max())
        plt.xlabel(r"$(4x^3 - 12Lx^2 + 9L^2x - L^3)/\mathrm{m}^3$")
        plt.ylabel(r"$D/\mathrm{mm}$")
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig("build/plot_{}_{}_{}_2.pdf".format(material, querschnitt, einspannung))
        plt.close()

        print("""
        ------------------------------------------------------------------------
        Material:                   {}
        Querschnitt:                {}
        Einspannung:                {}
        Durchmesser d:              {} ± {} mm
        Länge L:                    {} cm
        Masse M:                    {} kg
        Flächenträgheitsmoment I:   {} ± {} mm^4
        Elastizitätsmodul E1:  {} ± {} N/m^2
        Elastizitätsmodul E2:  {} ± {} N/m^2
        Elastizitätsmodul E:   {} ± {} N/m^2
        Steigung m1:           {} ± {}
        Steigung m2:           {} ± {}
        Ordinatenabschnitt b1: {} ± {}
        Ordinatenabschnitt b2: {} ± {}
        ------------------------------------------------------------------------
        """.format(material, querschnitt, einspannung, d*1e3, 0, L*1e2, M, I*1e12, 0, E1*1e0, ΔE1*1e0, E2*1e0, ΔE2*1e0, E*1e0, ΔE*1e0, m1, Δm1, m2, Δm2, b1, Δb1, b2, Δb2))

'''
############################################################################
# Test mit Messwerten von Philipp Leser
# Aluminium, kreisförmiger Querschnitt, einseitig eingespannt
# Daten einlesen
x, D = np.loadtxt("data/daten_aluminium_quadratisch_beidseitig.txt", unpack=True)
d = np.loadtxt("data/daten_aluminium_quadratisch_durchmesser.txt", unpack=True)
L = 55.3 #[cm]
M = 4.6944 #[kg]
# Auswertung
d *= 1e-3
L *= 1e-2
x *= 1e-2
D *= 1e-6
auswertung("Aluminium", "quadratisch", "beidseitig", x, D, d, L, M)
############################################################################
'''
# Aluminium, kreisförmiger Querschnitt, einseitig eingespannt
# Daten einlesen
x, D = np.loadtxt("Messing.txt", unpack=True)
d = 10 #[mm]
L = 40.70 #[cm]
M = 2.3606 #[kg]
# Auswertung
d *= 1e-3
L *= 1e-2
x *= 1e-2
D *= 1e-6
auswertung("Messing", "quadratisch", "einseitig", x, D, d, L, M)

# Messing, quadratischer Querschnitt, einseitig eingespannt
# Daten einlesen
x, D = np.loadtxt("alurund.txt", unpack=True)
d = 10 #[mm]
L = 34.8 #[cm]
M = 1.1926 #[kg]
# Auswertung
d *= 1e-3
L *= 1e-2
x *= 1e-2
D *= 1e-6
auswertung("Aluminium", "kreisfoermig", "einseitig", x, D, d, L, M)

# Stahl, quadratischer Querschnitt, beidseitig eingespannt
# Daten einlesen
x, D = np.loadtxt("alueckig.txt", unpack=True)
d = 10 #[mm]
L = 55.3 #[cm]
M = 0 #[kg]
# Auswertung
d *= 1e-3
L *= 1e-2
x *= 1e-2
D *= 1e-6
auswertung("Aluminium", "quadratisch", "beidseitig" , x, D, d, L, M)

x, D = np.loadtxt("alueckig2.txt", unpack=True)
d = 10 #[mm]
L = 55.3 #[cm]
M = 3.5312 #[kg]
# Auswertung
d *= 1e-3
L *= 1e-2
x *= 1e-2
D *= 1e-6
auswertung("Aluminium", "quadratisch", "beidseitig", x, D, d, L, M)
