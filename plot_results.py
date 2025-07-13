#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 11:14:26 2025

@author: kissami
"""

import numpy as np
from scipy.optimize import fsolve

def exact_SWE(x, h1, h2, g, t, xm):
    def racine(h2, h1, g):
        c2 = np.sqrt(h2 * g)
        c1 = np.sqrt(h1 * g)
        
        def eq(cm):
            return 8 * c1**2 * cm**2 * (c2 - cm)**2 - (cm**2 - c1**2)**2 * (cm**2 + c1**2)
        
        cm = fsolve(eq, c2)[0]
        hm = cm**2 / g
        return cm, hm

    c2 = np.sqrt(g * h2)
    cm, hm = racine(h2, h1, g)
    um = 2 * (c2 - cm)
    vc = hm * um / (hm - h1)
    xa = -t * c2
    xb = (2 * c2 - 3 * cm) * t
    xc = vc * t

    h = np.zeros_like(x)
    u = np.zeros_like(x)

    for k in range(len(x)):
        if x[k] - xm <= xa:
            h[k] = h2
            u[k] = 0
        elif xa < x[k] - xm <= xb:
            h[k] = 1 / (9 * g) * (2 * c2 - (x[k] - xm) / t)**2
            u[k] = 2 / 3 * (c2 + (x[k] - xm) / t)
        elif xb < x[k] - xm <= xc:
            h[k] = hm
            u[k] = um
        else:  # x[k] - xm > xc
            h[k] = h1
            u[k] = 0

    return h, u


# Paramètres du test
    

if __name__ == "__main__":
    import matplotlib.pyplot as plt


    h_left = 5.0  
    h_right = 2.0  
    g = 1      
    t = 3.0        
    xm = 0.0     

    # === Comparaison avec solution numérique (si disponible) ===
    data = np.loadtxt("output_mpi.dat")
    print(data)
    x_num = data[:, 0]
    h_num = data[:, 1]
    u_num = data[:, 2]
    
    # Domaine 
    x = np.linspace(-10, 10, x_num.shape[0])

    # Calcul de la solution
    h, u = exact_SWE(x, h1=h_right, h2=h_left, g=g, t=t, xm=xm)


    plt.figure(figsize=(10, 4))
    plt.plot(x, h, 'k--', linewidth=2, label='h exacte')
    plt.plot(x_num, h_num, 'b-', markersize=4, label='h numérique MPI')
    plt.plot(x, u, 'k:', linewidth=2, label='u exacte')
    plt.plot(x_num, u_num, 'r-', markersize=4, label='u numérique MPI')
    plt.xlabel("x")
    plt.ylabel("h(x), u(x)")
    plt.legend()
    plt.title(f"Dam Break: Numérique vs Exact à t={t:.1f} s")
    plt.grid()
    plt.tight_layout()

    plt.savefig("dam_break_comparison.png", dpi=300)
    # plt.show()

