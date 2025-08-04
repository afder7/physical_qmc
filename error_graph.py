from functools import reduce
from operator import mul
from math import *
import matplotlib.pyplot as plt
import random


n = 100
d = 5
f = lambda x: sqrt(18 / d) * (sum([sqrt(j) for j in x]) - 2 * d / 3)

with open(f"generation/modeled_num_{n}_dim_{d}.txt", "r") as fi:
    modeled_points = [[float(y) for y in x.split()] for x in fi.read().split("\n")[:-1]]

with open(f"generation/random_num_{n}_dim_{d}.txt", "r") as fi:
    random_points = [[float(y) for y in x.split()] for x in fi.read().split("\n")[:-1]]

gen = []
rnd = []
it_n = n
it = range(10, it_n + 1)
for i in range(10, it_n + 1):
    pt = [random.choice(modeled_points) for _ in range(i)]
    s = 0
    for x in pt:
        s += f(x)
    gen.append(log(abs(s / i), e))

    pt = [[random.random() for __ in range(d)] for _ in range(i)]
    s = 0
    for x in pt:
        s += f(x)
    rnd.append(log(abs(s / i), e))

    if i % 500 == 0:
        print(i)

print(gen[-1])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))  # 1 строка, 2 столбца

# Первый график (слева)
ax1.plot(it, gen, 'b-', label='Generated points', color='blue', lw=0.7)
ax1.set_title('Generated points e(N)')
ax1.set_xlabel('N, number of points')
ax1.set_ylabel('ln e, error of calculation')
# ax1.set_ylim(-0.3, 0.3)

# Второй график (справа)
ax2.plot(it, rnd, 'r-', label='Random points', color='green', lw=0.7)
ax2.set_title('Random points e(N)')
ax2.set_xlabel('N, number of points')
ax2.set_ylabel('ln e, error of calculation')
# ax2.set_ylim(-0.3, 0.3)

# Общие настройки
for ax in (ax1, ax2):
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()
