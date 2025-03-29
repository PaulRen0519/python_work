# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 14:00:16 2025

@author: admin
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
def laplacian(Z):
    Ztop = Z[0:-2, 1:-1]
    Zleft = Z[1:-1, 0:-2]
    Zbottom = Z[2:, 1:-1]
    Zright = Z[1:-1, 2:]
    Zcenter = Z[1:-1, 1:-1]
    return (Ztop + Zleft + Zbottom + Zright - 4*Zcenter)
Nx, Ny = 128, 128
Nt = 1000
f = 0.054
k = 0.062
Du = 1.0
Dv = 0.5
u = np.ones((Nx, Ny))
v = np.zeros((Nx, Ny))
r = 20


