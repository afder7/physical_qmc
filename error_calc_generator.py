import random
from math import *
from scipy.stats import qmc
import matplotlib.pyplot as plt
import signal
import copy
import numpy as np


number_of_points = 200
dimension = 100
iter_count = 1000
delta_t = 1 / number_of_points
kappa = 200

it = 0

G = 1
m = 1  # actual value is calculated when the points are scattered
c = 1  # actual value is calculated when the points are scattered
A0 = m / delta_t ** 2 + c / (2 * delta_t)
A1 = 2 * m / delta_t ** 2
A2 = m / delta_t ** 2 - c / (2 * delta_t)

cur_U = 0

counter = 0
x_prev = []
x_0 = []

potential_energy_list = []
l2_discrepancy_list = []
wd_discrepancy_list = []
star_discrepancy_list = []

f = [[]]
skip_iterations = 10

best_l2 = 10000
best_points = []


def dist(x_i, x_j):
    dim = len(x_i)
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


# points is a set of num points in a given dim
# U_{1,1}
def potential_energy(points):
    num = gnum
    dim = dimension

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


def force(points):
    num = gnum
    dim = dimension
    re = [[0 for __ in range(dim)] for _ in range(num)]

    for i in range(num):
        for j in range(num):
            if j != i:
                d = dist(points[i], points[j])
                for k in range(dim):
                    sign = -1 if points[i][k] < points[j][k] else 1
                    delta = abs(points[i][k] - points[j][k])
                    a_ijk = sign * delta * (1 - delta) * (1 - 2 * delta)
                    re[i][k] += -G * a_ijk / (d ** 3)

    return re


class MovingPoints:
    def __init__(self, point_count=number_of_points, dim=dimension):
        self.num = point_count
        self.dim = dim
        self.points = []
        self.velocities = []

        self.create_points()

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_U

        for i in range(self.num):
            point = [random.random() for _ in range(self.dim)]
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        U = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(U) * delta_t ** 2 / scalar
        c = -abs(U) * sqrt(kappa) * delta_t / scalar
        cur_U = U

        # calculating x_1
        global f
        f = force(self.points)
        for i in range(self.num):
            for k in range(self.dim):
                self.points[i][k] = self.points[i][k] + 0.5 * (1 / m) * f[i][k] * delta_t ** 2
                if self.points[i][k] > 1:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k])
                elif self.points[i][k] < 0:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k]) + 1

    def update(self):
        global x_prev, f
        if it % skip_iterations == 0:
            global best_l2, best_points
            f = force(self.points)
            l2_disc = qmc.discrepancy(self.points, method='L2-star')
            if l2_disc < best_l2:
                best_points = copy.deepcopy(self.points)
                best_l2 = l2_disc

        new_x_prev = copy.deepcopy(self.points)
        for i in range(self.num):
            for k in range(self.dim):
                x_j = self.points[i][k]
                self.points[i][k] = (-f[i][k] + A1 * x_j - A2 * x_prev[i][k]) / A0
                if self.points[i][k] > 1:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k])
                elif self.points[i][k] < 0:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k]) + 1

        x_prev = copy.deepcopy(new_x_prev)


func = lambda x: sqrt(18 / dimension) * (sum([sqrt(j) for j in x]) - 2 * dimension / 3)
gen = []

for gnum in range(2, number_of_points + 1):
    app = MovingPoints(dim=dimension, point_count=gnum)
    for _ in range(iter_count):
        # print("Iteration        ", it)
        app.update()
        it += 1
    it = 0

    s = 0
    for xi in best_points:
        s += func(xi)
    gen.append(log(abs(s / gnum), e))

    print(gnum)

x_ax = np.array(range(2, number_of_points + 1))
gen = np.array(gen)
k, b = np.polyfit(x_ax, gen, deg=1)  # deg=1 for linear fit
trendline = k * x_ax + b

plt.plot(x_ax, trendline, 'r--', label=f'Trendline: y = {k:.2f}x + {b:.2f}')
plt.plot(x_ax, gen, label='Generated points error trend', linewidth=1.2, color="blue")

plt.xlabel("N, number of points")
plt.ylabel("e, error of calculation")
plt.legend()
plt.grid(True)
plt.show()
