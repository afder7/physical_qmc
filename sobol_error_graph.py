from functools import reduce
from operator import mul
from math import *
import matplotlib.pyplot as plt
from scipy.stats import qmc
import numpy as np


n = 200
d = 100
f = lambda x: sqrt(18 / d) * (sum([sqrt(j) for j in x]) - 2 * d / 3)

sobol_engine = qmc.Sobol(d=d)  # scramble=True для рандомизации
sobol_scramble_engine = qmc.Sobol(d=d, scramble=True)

gen = []
rnd = []
it_n = n
it = range(2, it_n + 1)
for i in range(2, it_n + 1):
    pt = list(sobol_engine.random(i))
    s = 0
    for x in pt:
        s += f(x)
    gen.append(log(abs(s / i), e))

    if i % 500 == 0:
        print(i)

x_ax = np.array(range(2, n + 1))
gen = np.array(gen)
k, b = np.polyfit(x_ax, gen, deg=1)  # deg=1 for linear fit
trendline = k * x_ax + b

plt.plot(it, trendline, 'r--', label=f'Trendline: y = {k:.2f}x + {b:.2f}')
plt.plot(it, gen, 'b-', label='Sobol points', color='blue', lw=0.7)
plt.title('Sobol points ln e(N)')
plt.xlabel('N, number of points')
plt.ylabel('ln e, error of calculation')
# ax1.set_ylim(-0.2, 0.2)

plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
