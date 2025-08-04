import random
from math import *
from scipy.stats import qmc
import matplotlib.pyplot as plt
import signal
import copy
import numpy as np


number_of_points = 256
dimension = 1000
iter_count = 500
skip_iterations = 10
delta_t = 1 / number_of_points
kappa = 200


eps = 1
m = 1  # actual value is calculated when the points are scattered
c = 1  # actual value is calculated when the points are scattered
A0 = m / delta_t ** 2 + c / (2 * delta_t)
A1 = 2 * m / delta_t ** 2
A2 = m / delta_t ** 2 - c / (2 * delta_t)

cur_u = 0
counter = 0
x_prev = []
x_0 = []

l2_discrepancy_list = []

f = [[]]
it = 0

random_l2_disc = 10000
best_l2_disc = 10000
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
    num = number_of_points
    dim = dimension

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


def force(points):
    num = number_of_points
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
                    re[i][k] += -eps * a_ijk / (d ** 3)

    return re


class PointGenerator:
    def __init__(self):
        self.num = number_of_points
        self.dim = dimension
        self.points = []
        self.velocities = []

        self.create_points()

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_u, l2_discrepancy_list

        for i in range(self.num):
            point = [random.random() for _ in range(self.dim)]
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        global random_l2_disc, best_l2_disc, best_points
        save_points(self.points, "random")
        random_l2_disc = qmc.discrepancy(self.points, method='L2-star')
        l2_discrepancy_list.append(random_l2_disc)
        if random_l2_disc < best_l2_disc:
            best_points = copy.deepcopy(self.points)
            best_l2_disc = random_l2_disc

        u = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(u) * delta_t ** 2 / scalar
        c = -abs(u) * sqrt(kappa) * delta_t / scalar
        cur_u = u

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
            f = force(self.points)
        else:
            f = [[y / skip_iterations for y in x] for x in f]

        new_x_prev = self.points
        for i in range(self.num):
            for k in range(self.dim):
                x_j = self.points[i][k]
                self.points[i][k] = (-f[i][k] + A1 * x_j - A2 * x_prev[i][k]) / A0
                if self.points[i][k] > 1:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k])
                elif self.points[i][k] < 0:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k]) + 1
        x_prev = new_x_prev.copy()

        if it % skip_iterations == 0:
            global best_l2_disc, best_points, l2_discrepancy_list
            l2_disc = qmc.discrepancy(self.points, method='L2-star')
            l2_discrepancy_list.append(l2_disc)
            if l2_disc < best_l2_disc:
                best_points = copy.deepcopy(self.points)
                best_l2_disc = l2_disc


def save_points(points, save_type):
    with open(f"generation/{save_type}_num_{number_of_points}_dim_{dimension}.txt", "w") as file:
        for pt in points:
            for cd in pt:
                if cd > 1:
                    file.write(str(cd - int(cd)) + " ")
                elif cd < 0:
                    file.write(str(cd - int(cd) + 1) + " ")
                else:
                    file.write(str(cd) + " ")
            file.write("\n")


app = PointGenerator()


def build_graph():
    plt.figure(figsize=(12, 7))
    plt.title(f"num={number_of_points}, dim={dimension}, it={it}, skip={skip_iterations}")
    plt.plot([skip_iterations * x for x in range(len(l2_discrepancy_list))], l2_discrepancy_list, 'm--', label='L2')
    plt.legend()
    plt.tight_layout()
    plt.show()


def save_on_exit():
    print(f"\nSaving results of current iteration #{it}")
    save_points(app.points, "modeled")
    save_points(best_points, "best")

    sobol_engine = qmc.Sobol(d=dimension)  # scramble=True для рандомизации
    sobol_points = np.array(sobol_engine.random(number_of_points))
    modeled_l2_disc = qmc.discrepancy(app.points, method='L2-star')
    print("Random L2  ", random_l2_disc)
    print("Modeled L2 ", modeled_l2_disc)
    print("Best L2    ", best_l2_disc)
    print("Sobol L2   ", qmc.discrepancy(sobol_points, method='L2-star'))
    build_graph()

    exit(0)


signal.signal(signal.SIGINT, save_on_exit)
signal.signal(signal.SIGTERM, save_on_exit)


for _ in range(iter_count):
    if it % 20 == 0:
        print("Iteration        ", it)
    app.update()
    it += 1

save_on_exit()
