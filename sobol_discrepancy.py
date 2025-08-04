import numpy as np
from scipy.stats import qmc
from math import *
import matplotlib.pyplot as plt
import tkinter as tk


num = 4096
dim = 2


def dist(x_i, x_j):
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


def potential_energy(points):
    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


sobol_engine = qmc.Sobol(d=dim)  # scramble=True для рандомизации
points = np.array(sobol_engine.random(num))
print("Sobol L2-discrepancy ", qmc.discrepancy(points, method='L2-star'))
print("Sobol Potential Energy", potential_energy(points))
