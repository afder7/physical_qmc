import numpy as np
from scipy.stats import qmc

num = 100
dim = 2

sobol_engine = qmc.Sobol(d=dim)  # scramble=True для рандомизации
points = np.array(sobol_engine.random(num))
print("Sobol L2-discrepancy     ", format(qmc.discrepancy(points, method='L2-star'), '.8f'))
print("Generated L2-discrepancy ", format(0.009894785122310719, '.8f'))
