import numpy as np
from scipy.stats import qmc

num = 128
dim = 20

sobol_engine = qmc.Sobol(d=dim)  # scramble=True для рандомизации
points = np.array(sobol_engine.random(num))
print("Sobol L2-discrepancy ", format(qmc.discrepancy(points, method='L2-star'), '.8f'))
print("Sobol L2-discrepancy ", format(3.944450747283431e-05, '.8f'))
